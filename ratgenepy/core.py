"""Main RatGene entry point."""

from __future__ import annotations

import numpy as np

from ._utils import _as_str, _first_match
from .exchange import introExchange
from .integrate import uniteGeneModel
from .model import Model
from .ratio import ratioGene, ratioMethod
from .reduction import rmRedundancy
from .solvers import optimizeCbModel

def RatGene(model, targetMet, biomass=None, carbon="EX_glc__D_e",
            oxygen="EX_o2_e", LBbiomass=0.05, LBcarbon=-15.0,
            LBoxygen=-15.0, maxLoop=1000, timeLimit=float("inf"), pool=10,
            type=1, size=1, addition=None, backend="scip"):
    """Ratio-based gene/reaction modification strategy search.

    Parameters follow ``RatGene.m``.  ``type`` is 1 for genes and 0 for
    reactions; ``size`` is 1 when the reduction procedure is applied.
    ``backend`` selects the LP/MILP solver: ``"scip"`` (default),
    ``"cplex"`` or ``"gurobi"``.

    Returns
    -------
    (x_target, knockouts) or (x_target, knockouts, additions)
    """
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

    if type == 1:
        if addition_problem:
            model1 = uniteGeneModel(model, addition)
        else:
            model1 = model.copy()

        model1.lb[id_biomass] = LBbiomass
        model1.lb[id_carbon] = LBcarbon
        model1.lb[id_oxygen] = LBoxygen
        model1.ub[id_carbon] = 0
        model1.ub[id_oxygen] = 0

        new_model, id_target, TMPR = introExchange(
            model1, id_biomass, [id_carbon, id_oxygen], metID,
            backend=backend
        )
        if TMPR <= 0:
            raise ValueError("theoretical maximum production rate is not positive")

        opt = optimizeCbModel(model1, backend=backend)
        TMGR = opt.f

        x_target, _, knockout = ratioGene(
            new_model, id_biomass, id_target, TMGR, maxLoop,
            TMPR / (maxLoop * LBbiomass), pool, timeLimit,
            backend=backend
        )

        if np.asarray(knockout).size == 0:
            return (x_target, [], []) if addition_problem else (x_target, [])

        if size == 0:
            newKnockout = knockout
        else:
            split = ([len(model.rxns), len(model.genes)]
                     if addition_problem else [0, 0])
            _, _, newKnockout = rmRedundancy(
                new_model, knockout, id_target, id_biomass, split,
                backend=backend
            )

        if not addition_problem:
            return (x_target,
                    [new_model.genes[i] for i in np.where(np.asarray(newKnockout) == 0)[0]])

        n_core_genes = len(model.genes)
        delete = [model.genes[i] for i in np.where(
            np.asarray(newKnockout)[:n_core_genes] == 0)[0]]
        addGene = list(new_model.genes[n_core_genes:])
        add_values = np.asarray(newKnockout)[n_core_genes:]
        additions = [addGene[i] for i in np.where(add_values == 1)[0]]
        return x_target, delete, additions

    # reaction deletion (or reaction deletion/addition)
    if addition_problem:
        model1 = uniteGeneModel(model, addition)
    else:
        model1 = model.copy()

    model1.lb[id_biomass] = LBbiomass
    model1.lb[id_carbon] = LBcarbon
    model1.lb[id_oxygen] = LBoxygen
    model1.ub[id_carbon] = 0
    model1.ub[id_oxygen] = 0

    new_model, id_target, TMPR = introExchange(
        model1, id_biomass, [id_carbon, id_oxygen], metID,
        backend=backend
    )
    if TMPR <= 0:
        raise ValueError("theoretical maximum production rate is not positive")

    x_target, _, newKnockout = ratioMethod(
        new_model, id_biomass, id_target, maxLoop,
        TMPR / (maxLoop * LBbiomass), timeLimit,
        backend=backend
    )
    newKnockout = np.asarray(newKnockout)

    if not addition_problem:
        return x_target, [new_model.rxns[i] for i in np.where(newKnockout == 0)[0]]

    n_core_rxns = len(model.rxns)
    delete = [model.rxns[i] for i in np.where(newKnockout[:n_core_rxns] == 0)[0]]
    if id_target == len(new_model.rxns) - 1:
        addRxn = list(new_model.rxns[n_core_rxns:len(new_model.rxns) - 1])
        add_values = newKnockout[n_core_rxns:len(new_model.rxns) - 1]
    else:
        addRxn = list(new_model.rxns[n_core_rxns:])
        add_values = newKnockout[n_core_rxns:]
    additions = [addRxn[i] for i in np.where(add_values == 1)[0]]
    return x_target, delete, additions
