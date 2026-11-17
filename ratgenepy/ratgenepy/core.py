"""Main RatGene entry point."""

from __future__ import annotations

import numpy as np

from ._utils import _as_str, _first_match
from .changelog import ChangeLog, gpr_genes, write_change_log
from .exchange import _normalize_target_type, introExchange
from .gpr import verifyRatioGene
from .integrate import uniteGeneModel
from .model import Model
from .ratio import ratioGene, ratioMethod
from .reduction import rmRedundancy
from .solvers import optimizeCbModel

#: Change log of the most recent :func:`RatGene` call.  This is provided as a
#: convenient way to inspect the log when the default return value is kept.
LAST_CHANGE_LOG = None


class RatGeneResult(tuple):
    """Tuple returned by :func:`RatGene`.

    It behaves exactly like the historical result tuple (and therefore keeps
    unpacking such as ``x_target, knockouts = RatGene(...)`` working) while
    also exposing the generated change log as ``result.change_log``.
    """

    def __new__(cls, values, change_log: ChangeLog):
        result = super().__new__(cls, tuple(values))
        result.change_log = change_log
        return result

    @property
    def log(self):
        """Alias for :attr:`change_log`."""
        return self.change_log

    @property
    def target_info(self):
        """Metadata about the selected target outlet, if available.

        The dictionary identifies whether an ``EX_`` secretion outlet or a
        ``DM_`` intracellular demand was used (keys ``kind``, ``role``,
        ``reaction_id``, ``introduced``, ``opened``, ``ambiguous`` ...).
        """
        return self.change_log.metadata.get("target", {})


def _set_bounds(model, index: int, lb: float, ub: float, log: ChangeLog,
                reason: str):
    """Set bounds and record the change when it is not a no-op."""
    index = int(index)
    old = (float(model.lb[index]), float(model.ub[index]))
    model.lb[index] = float(lb)
    model.ub[index] = float(ub)
    log.add_bounds_change(model, index, old, (lb, ub), reason)


def _log_gene_deletions(log: ChangeLog, model: Model, gene_vector,
                        reaction_limit: int, *, label: str = "final strategy"):
    """Record reactions disabled by a gene strategy."""
    gene_vector = np.asarray(gene_vector).ravel()
    reaction_status = verifyRatioGene(model, gene_vector)
    deleted_genes = [
        model.genes[i] for i in range(min(len(gene_vector), len(model.genes)))
        if gene_vector[i] == 0
    ]
    for index in np.where(reaction_status == 1)[0]:
        index = int(index)
        if index >= reaction_limit:
            continue
        old = (float(model.lb[index]), float(model.ub[index]))
        if old == (0.0, 0.0):
            continue
        rule_genes = gpr_genes(model.grRules[index]) if index < len(model.grRules) else []
        genes = [gene for gene in deleted_genes if gene in rule_genes]
        if not genes:
            genes = rule_genes
        if genes:
            reason = f"gene deletion in {label}: {', '.join(genes)}"
        else:
            reason = f"gene deletion in {label}"
        log.add_bounds_change(
            model, index, old, (0.0, 0.0), reason,
            change_type="reaction_knockout", role="gene_deletion",
            genes=genes,
        )


def _log_gene_additions(log: ChangeLog, model: Model, gene_vector,
                        n_core_genes: int, n_core_rxns: int):
    """Record reactions enabled by genes taken from an addition model."""
    gene_vector = np.asarray(gene_vector).ravel()
    if len(gene_vector) <= n_core_genes:
        return
    added_indices = [
        i for i in range(n_core_genes, min(len(gene_vector), len(model.genes)))
        if gene_vector[i] == 1
    ]
    if not added_indices:
        return
    added_genes = [model.genes[i] for i in added_indices]

    # Compare the GPR state with all core genes present and all added genes
    # absent against the final strategy.  Only reactions that become active
    # because of an added gene are reported.
    before = np.ones(len(model.genes))
    before[n_core_genes:] = 0.0
    status_before = verifyRatioGene(model, before)
    status_after = verifyRatioGene(model, gene_vector)

    for index in np.where((status_before == 1) & (status_after == 0))[0]:
        index = int(index)
        rule_genes = gpr_genes(model.grRules[index]) if index < len(model.grRules) else []
        genes = [gene for gene in added_genes if gene in rule_genes]
        if not genes:
            continue
        new = (float(model.lb[index]), float(model.ub[index]))
        if new == (0.0, 0.0):
            continue
        introduced = index >= n_core_rxns
        old = (0.0, 0.0) if introduced else new
        rule = model.grRules[index] if index < len(model.grRules) else ""
        log.add_bounds_change(
            model, index, old, new,
            f"gene addition enables reaction (GPR: {rule})",
            change_type="gene_addition", role="gene_addition",
            genes=genes, introduced=introduced, force=True,
        )


def _log_reaction_deletions(log: ChangeLog, model: Model, knockout,
                            reaction_limit: int):
    knockout = np.asarray(knockout).ravel()
    for index in np.where(knockout == 0)[0]:
        index = int(index)
        if index >= reaction_limit:
            continue
        old = (float(model.lb[index]), float(model.ub[index]))
        if old == (0.0, 0.0):
            continue
        log.add_bounds_change(
            model, index, old, (0.0, 0.0),
            "reaction deletion in final strategy",
            change_type="reaction_knockout", role="reaction_deletion",
        )


def _log_reaction_additions(log: ChangeLog, model: Model, knockout,
                            n_core_rxns: int, id_target: int):
    knockout = np.asarray(knockout).ravel()
    for index in range(n_core_rxns, len(knockout)):
        if index == id_target or knockout[index] != 1:
            continue
        new = (float(model.lb[index]), float(model.ub[index]))
        if new == (0.0, 0.0):
            continue
        log.add_bounds_change(
            model, index, (0.0, 0.0), new,
            "reaction addition in final strategy",
            change_type="reaction_addition", role="reaction_addition",
            introduced=True, force=True,
        )


def _finalize(result, log: ChangeLog, return_changelog: bool,
              changelog, change_log_file):
    global LAST_CHANGE_LOG
    LAST_CHANGE_LOG = log
    RatGene.last_change_log = log

    if changelog is not None and changelog is not log:
        try:
            changelog.extend(log)
        except AttributeError as exc:  # pragma: no cover - user error
            raise TypeError("changelog must be a ChangeLog or a list") from exc

    if change_log_file is not None:
        write_change_log(log, change_log_file)

    result = tuple(result)
    if return_changelog:
        result = result + (log,)
    return RatGeneResult(result, log)


def RatGene(model, targetMet, biomass=None, carbon="EX_glc__D_e",
            oxygen="EX_o2_e", LBbiomass=0.05, LBcarbon=-15.0,
            LBoxygen=-15.0, maxLoop=1000, timeLimit=float("inf"), pool=10,
            type=1, size=1, addition=None, backend="scip",
            integerTolerance=None, target_type="auto",
            return_changelog=False, changelog=None, change_log_file=None,
            **kwargs):
    """Ratio-based gene/reaction modification strategy search.

    Parameters follow ``RatGene.m``.  ``type`` is 1 for genes and 0 for
    reactions; ``size`` is 1 when the reduction procedure is applied.
    ``backend`` selects the LP/MILP solver: ``"scip"`` (default),
    ``"cplex"`` or ``"gurobi"``.  ``integerTolerance`` sets the MILP
    integrality tolerance; when omitted, each backend uses its smallest
    supported value.

    ``target_type`` selects the target production outlet explicitly:
    ``"auto"`` (prefer a secretion outlet ``EX_``, then an intracellular
    demand ``DM_``), ``"secretion"`` or ``"demand"``.

    A change log is always generated.  Pass ``return_changelog=True`` to
    append it to the returned tuple, ``changelog=[]`` to fill an external
    list, or ``change_log_file="changes.json"`` (``.csv`` and text are also
    supported) to write it to disk.  The most recent log is also available as
    ``RatGene.last_change_log`` and :data:`ratgenepy.core.LAST_CHANGE_LOG`.
    The selected target outlet (``EX_`` secretion or ``DM_`` demand) is
    recorded in ``change_log.metadata["target"]`` and exposed as
    ``result.target_info``.

    Returns
    -------
    (x_target, knockouts) or (x_target, knockouts, additions)
    (+ change_log when ``return_changelog`` is true)
    """
    # Accept a couple of descriptive aliases without changing the public API.
    if "return_change_log" in kwargs:
        return_changelog = kwargs.pop("return_change_log")
    if "return_log" in kwargs:
        return_changelog = kwargs.pop("return_log")
    if "returnChangeLog" in kwargs:
        return_changelog = kwargs.pop("returnChangeLog")
    if "targetType" in kwargs:
        target_type = kwargs.pop("targetType")
    if "changeLogFile" in kwargs:
        if change_log_file is not None:
            raise TypeError("use either 'change_log_file' or 'changeLogFile', not both")
        change_log_file = kwargs.pop("changeLogFile")
    if "change_log" in kwargs:
        if changelog is not None:
            raise TypeError("use either 'changelog' or 'change_log', not both")
        value = kwargs.pop("change_log")
        if isinstance(value, bool):
            return_changelog = value
        else:
            changelog = value
    if "changeLog" in kwargs:
        if changelog is not None:
            raise TypeError("use either 'changelog' or 'changeLog', not both")
        value = kwargs.pop("changeLog")
        if isinstance(value, bool):
            return_changelog = value
        else:
            changelog = value
    if kwargs:
        unexpected = next(iter(kwargs))
        raise TypeError(f"RatGene() got an unexpected keyword argument {unexpected!r}")

    target_type = _normalize_target_type(target_type)
    return_changelog = bool(return_changelog)

    if isinstance(changelog, ChangeLog):
        log = changelog
        external_log = None
    else:
        log = ChangeLog()
        external_log = changelog

    if not isinstance(model, Model):
        model = Model.from_mat(model)

    # target metabolite
    if isinstance(targetMet, (list, tuple, np.ndarray)):
        targetMet = _as_str(np.asarray(targetMet).ravel()[0])
    targetMet = _as_str(targetMet)
    metID = _first_match(model.mets, targetMet)

    # biomass
    if biomass is not None:
        if isinstance(biomass, (list, tuple, np.ndarray)):
            biomass = _as_str(np.asarray(biomass).ravel()[0])
        id_biomass = _first_match(model.rxns, _as_str(biomass))
    else:
        nonEmpty = [i for i, name in enumerate(model.rxns) if "BIOMASS" in name]
        if len(nonEmpty) > 1:
            potential = [model.rxns[i] for i in nonEmpty]
            ni = [i for i, name in enumerate(potential) if "core" in name]
            id_biomass = nonEmpty[ni[0]] if ni else nonEmpty[0]
        elif len(nonEmpty) == 1:
            id_biomass = nonEmpty[0]
        else:
            raise ValueError("no biomass rxn in the model")

        c_index = np.where(model.c != 0)[0]
        if c_index.size:
            id_biomass = int(c_index[0])

    # carbon and oxygen sources
    id_carbon = _first_match(model.rxns, _as_str(carbon))
    id_oxygen = _first_match(model.rxns, _as_str(oxygen))

    if LBbiomass <= 0:
        raise ValueError("biomass lower bound should be positive")
    if type not in (0, 1):
        raise ValueError("type can only be 1-gene or 0-reaction")
    if size not in (0, 1):
        raise ValueError("size can only be 1 or 0")

    addition_problem = addition is not None
    if addition_problem and not isinstance(addition, Model):
        addition = Model.from_mat(addition)

    n_core_rxns = len(model.rxns)
    n_core_genes = len(model.genes)

    if type == 1:
        if addition_problem:
            model1 = uniteGeneModel(model, addition)
        else:
            model1 = model.copy()

        _set_bounds(
            model1, id_biomass, LBbiomass, float(model1.ub[id_biomass]), log,
            f"biomass lower bound set to LBbiomass={LBbiomass:g}",
        )
        _set_bounds(
            model1, id_carbon, LBcarbon, 0.0, log,
            f"carbon source uptake constraint (LBcarbon={LBcarbon:g}, ub=0)",
        )
        _set_bounds(
            model1, id_oxygen, LBoxygen, 0.0, log,
            f"oxygen source uptake constraint (LBoxygen={LBoxygen:g}, ub=0)",
        )

        new_model, id_target, TMPR, target_info = introExchange(
            model1, id_biomass, [id_carbon, id_oxygen], metID,
            backend=backend, target_type=target_type, return_info=True,
            change_log=log,
        )
        log.metadata["target"] = dict(target_info)
        if TMPR <= 0:
            raise ValueError("theoretical maximum production rate is not positive")

        opt = optimizeCbModel(model1, backend=backend)
        TMGR = opt.f

        x_target, _, knockout = ratioGene(
            new_model, id_biomass, id_target, TMGR, maxLoop,
            TMPR / (maxLoop * LBbiomass), pool, timeLimit,
            backend=backend, integerTolerance=integerTolerance
        )

        if np.asarray(knockout).size == 0:
            result = (x_target, [], []) if addition_problem else (x_target, [])
            return _finalize(result, log, return_changelog, external_log,
                             change_log_file)

        if size == 0:
            newKnockout = knockout
        else:
            split = ([n_core_rxns, n_core_genes]
                     if addition_problem else [0, 0])
            _, _, newKnockout = rmRedundancy(
                new_model, knockout, id_target, id_biomass, split,
                backend=backend, integerTolerance=integerTolerance
            )

        # Log the reaction-level consequences of the final gene strategy.
        _log_gene_deletions(log, new_model, newKnockout, n_core_rxns)
        if addition_problem:
            _log_gene_additions(log, new_model, newKnockout,
                                n_core_genes, n_core_rxns)

        if not addition_problem:
            result = (
                x_target,
                [new_model.genes[i]
                 for i in np.where(np.asarray(newKnockout) == 0)[0]],
            )
            return _finalize(result, log, return_changelog, external_log,
                             change_log_file)

        delete = [model.genes[i] for i in np.where(
            np.asarray(newKnockout)[:n_core_genes] == 0)[0]]
        addGene = list(new_model.genes[n_core_genes:])
        add_values = np.asarray(newKnockout)[n_core_genes:]
        additions = [addGene[i] for i in np.where(add_values == 1)[0]]
        result = (x_target, delete, additions)
        return _finalize(result, log, return_changelog, external_log,
                         change_log_file)

    # reaction deletion (or reaction deletion/addition)
    if addition_problem:
        model1 = uniteGeneModel(model, addition)
    else:
        model1 = model.copy()

    _set_bounds(
        model1, id_biomass, LBbiomass, float(model1.ub[id_biomass]), log,
        f"biomass lower bound set to LBbiomass={LBbiomass:g}",
    )
    _set_bounds(
        model1, id_carbon, LBcarbon, 0.0, log,
        f"carbon source uptake constraint (LBcarbon={LBcarbon:g}, ub=0)",
    )
    _set_bounds(
        model1, id_oxygen, LBoxygen, 0.0, log,
        f"oxygen source uptake constraint (LBoxygen={LBoxygen:g}, ub=0)",
    )

    new_model, id_target, TMPR, target_info = introExchange(
        model1, id_biomass, [id_carbon, id_oxygen], metID,
        backend=backend, target_type=target_type, return_info=True,
        change_log=log,
    )
    log.metadata["target"] = dict(target_info)
    if TMPR <= 0:
        raise ValueError("theoretical maximum production rate is not positive")

    x_target, _, newKnockout = ratioMethod(
        new_model, id_biomass, id_target, maxLoop,
        TMPR / (maxLoop * LBbiomass), timeLimit,
        backend=backend
    )
    newKnockout = np.asarray(newKnockout)

    if newKnockout.size:
        _log_reaction_deletions(
            log, new_model, newKnockout,
            n_core_rxns if addition_problem else len(new_model.rxns),
        )

    if not addition_problem:
        result = (x_target,
                  [new_model.rxns[i] for i in np.where(newKnockout == 0)[0]])
        return _finalize(result, log, return_changelog, external_log,
                         change_log_file)

    if newKnockout.size:
        _log_reaction_additions(log, new_model, newKnockout,
                                n_core_rxns, id_target)

    delete = [model.rxns[i] for i in np.where(newKnockout[:n_core_rxns] == 0)[0]]
    if id_target == len(new_model.rxns) - 1:
        addRxn = list(new_model.rxns[n_core_rxns:len(new_model.rxns) - 1])
        add_values = newKnockout[n_core_rxns:len(new_model.rxns) - 1]
    else:
        addRxn = list(new_model.rxns[n_core_rxns:])
        add_values = newKnockout[n_core_rxns:]
    additions = [addRxn[i] for i in np.where(add_values == 1)[0]]
    result = (x_target, delete, additions)
    return _finalize(result, log, return_changelog, external_log,
                     change_log_file)
