"""Detection or introduction of the target metabolite exchange reaction.

This module distinguishes two different kinds of target production outlets:

``secretion``
    An exchange reaction (``EX_...``) that exports the target metabolite.
    A virtual secretion outlet (``EX_target``) is created when no suitable
    exchange reaction exists.

``demand``
    An intracellular demand reaction (``DM_...``).  In a steady-state model
    the flux through a demand reaction is the rate at which the target is
    produced and consumed internally.

The distinction is explicit: :func:`introExchange` can return metadata that
identifies which kind of outlet was selected, and it never silently replaces
an exchange reaction with a demand reaction (or vice versa).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .changelog import ChangeLog, reaction_compartment
from .model import Model
from .solvers import _maximize

TARGET_TYPES = ("auto", "secretion", "demand")

_TARGET_TYPE_ALIASES = {
    "auto": "auto",
    "automatic": "auto",
    "secretion": "secretion",
    "secretion_outlet": "secretion",
    "secretionoutlet": "secretion",
    "exchange": "secretion",
    "ex": "secretion",
    "outlet": "secretion",
    "demand": "demand",
    "intracellular_demand": "demand",
    "intracellulardemand": "demand",
    "dm": "demand",
}


def _normalize_target_type(target_type) -> str:
    value = "auto" if target_type is None else str(target_type)
    value = value.strip().lower().replace("-", "_").replace(" ", "_")
    try:
        return _TARGET_TYPE_ALIASES[value]
    except KeyError as exc:
        raise ValueError(
            "target_type must be 'auto', 'secretion' (EX_ reactions) or "
            "'demand' (DM_ reactions)"
        ) from exc


def _target_info(model: Model, id_target: int, kind: str, *, introduced: bool,
                 ambiguous: bool, opened: bool, candidates: dict,
                 old_bounds, new_bounds) -> dict:
    """Build the metadata dictionary returned by :func:`introExchange`."""
    reaction_id = str(model.rxns[id_target])
    compartment = reaction_compartment(model, id_target)
    role = "secretion_outlet" if kind == "secretion" else "intracellular_demand"
    return {
        "target_id": int(id_target),
        "reaction_id": reaction_id,
        "kind": kind,
        "type": kind,
        "target_type": kind,
        "role": role,
        "introduced": bool(introduced),
        "ambiguous": bool(ambiguous),
        "opened": bool(opened),
        "candidates": dict(candidates),
        "compartment": compartment,
        "old_bounds": tuple(old_bounds),
        "new_bounds": tuple(new_bounds),
    }


def _is_target_outlet(model: Model, id_met: int, index: int) -> bool:
    """Return true when reaction ``index`` is a pure target outlet.

    A production outlet consumes only the target metabolite, so positive flux
    through it represents secretion (or an intracellular demand).
    """
    try:
        column = model.S[:, index]
    except Exception:  # pragma: no cover - malformed model
        return False
    if hasattr(column, "toarray"):
        column = column.toarray()
    column = np.asarray(column).ravel()
    if float(column[id_met]) >= -1e-12:
        return False
    # Exchange and demand outlets consume only their target metabolite.
    return int(np.count_nonzero(np.abs(column) > 1e-12)) == 1


def _find_outlet(model: Model, id_met: int, prefix: str,
                 exact: Optional[str] = None) -> Optional[int]:
    """Find an existing ``EX_``/``DM_`` outlet for the target metabolite."""
    names = [str(name) for name in model.rxns]
    if exact is not None and exact in names:
        index = names.index(exact)
        if _is_target_outlet(model, id_met, index):
            return index
    matches = [
        index for index, name in enumerate(names)
        if name.startswith(prefix) and _is_target_outlet(model, id_met, index)
    ]
    if not matches:
        return None
    # Prefer an outlet that is already open for positive production flux.
    for index in matches:
        if float(model.ub[index]) > 0:
            return index
    return matches[0]


def _open_outlet(model: Model, index: int, change_log: Optional[ChangeLog],
                 kind: str) -> bool:
    """Allow positive flux through an existing but blocked outlet."""
    if float(model.ub[index]) > 0:
        return False
    old_bounds = (float(model.lb[index]), float(model.ub[index]))
    new_bounds = (old_bounds[0], max(1000.0, old_bounds[0]))
    model.ub[index] = new_bounds[1]
    if change_log is not None:
        role = "secretion_outlet" if kind == "secretion" else "intracellular_demand"
        change_log.add_bounds_change(
            model, index, old_bounds, new_bounds,
            f"opened blocked target {kind} outlet so it can carry "
            f"positive production flux",
            change_type="bounds", role=role,
        )
    return True


def introExchange(input_model: Model, id_biomass: int, id_input: Sequence[int],
                  id_met: int, backend: str = "scip",
                  target_type: str = "auto", return_info: bool = False,
                  change_log: Optional[ChangeLog] = None):
    """Detect or add a target production outlet.

    Parameters
    ----------
    input_model, id_biomass, id_input, id_met, backend
        Unchanged from the original RatGene interface.
    target_type
        ``"auto"`` (default) prefers an existing secretion outlet (``EX_``),
        then an existing intracellular demand (``DM_``), and finally creates a
        virtual secretion outlet.  ``"secretion"`` forces a secretion outlet
        (creating ``EX_target`` if necessary) and ``"demand"`` forces an
        intracellular demand (creating ``DM_target`` if necessary).
    return_info
        When true, a fourth value is returned with metadata about the selected
        outlet (reaction id, kind/role, compartment, bounds, whether it was
        introduced or opened, and candidates).
    change_log
        Optional :class:`ratgenepy.changelog.ChangeLog` to which a newly
        introduced target reaction -- or the bounds change used to reopen a
        blocked outlet -- is recorded.

    Returns
    -------
    ``(model, id_target, TMPR)`` or
    ``(model, id_target, TMPR, info)`` when ``return_info`` is true.
    """
    model = input_model.copy()
    target_type = _normalize_target_type(target_type)
    id_input = list(np.asarray(id_input).ravel())

    if len(id_input) == 1:
        id_carbon = int(id_input[0])
        id_oxygen = -1
    elif len(id_input) == 2:
        id_carbon = int(id_input[0])
        id_oxygen = int(id_input[1])
    else:
        raise ValueError("Input id error. At most 2 allowed.")

    name_met = str(model.mets[id_met])
    ex_rxn = "EX_" + name_met
    dx_rxn = "DM_" + name_met

    ex_id = _find_outlet(model, id_met, "EX_", exact=ex_rxn)
    dm_id = _find_outlet(model, id_met, "DM_", exact=dx_rxn)
    candidates = {"secretion": ex_id, "demand": dm_id}

    introduced = False
    ambiguous = False
    opened = False

    if target_type == "auto":
        if ex_id is not None:
            id_target = ex_id
            kind = "secretion"
            ambiguous = dm_id is not None
        elif dm_id is not None:
            id_target = dm_id
            kind = "demand"
        else:
            id_target = _add_target_reaction(model, id_met, "secretion")
            kind = "secretion"
            introduced = True
    elif target_type == "secretion":
        if ex_id is not None:
            id_target = ex_id
            kind = "secretion"
            ambiguous = dm_id is not None
        else:
            # The caller explicitly requested a secretion outlet; create one
            # even if an intracellular demand reaction exists.
            id_target = _add_target_reaction(model, id_met, "secretion")
            kind = "secretion"
            introduced = True
            ambiguous = dm_id is not None
    else:  # demand
        if dm_id is not None:
            id_target = dm_id
            kind = "demand"
            ambiguous = ex_id is not None
        else:
            # The caller explicitly requested an intracellular demand; create
            # one even if a secretion outlet exists.
            id_target = _add_target_reaction(model, id_met, "demand")
            kind = "demand"
            introduced = True
            ambiguous = ex_id is not None

    if id_target == id_carbon or id_target == id_oxygen:
        raise ValueError("target metabolite cannot be the carbon or oxygen source")

    if introduced:
        old_bounds = (0.0, 0.0)
        new_bounds = (float(model.lb[id_target]), float(model.ub[id_target]))
        if change_log is not None:
            role = "secretion_outlet" if kind == "secretion" else "intracellular_demand"
            if kind == "secretion":
                reason = (
                    f"introduced a virtual secretion outlet for target "
                    f"metabolite {name_met!r}"
                )
            else:
                reason = (
                    f"introduced a virtual intracellular demand for target "
                    f"metabolite {name_met!r}"
                )
            change_log.add_introduced_reaction(
                model, id_target, reason, role=role, old_bounds=old_bounds
            )
    else:
        old_bounds = (float(model.lb[id_target]), float(model.ub[id_target]))
        opened = _open_outlet(model, id_target, change_log, kind)
        new_bounds = (float(model.lb[id_target]), float(model.ub[id_target]))

    info = _target_info(
        model, id_target, kind, introduced=introduced, ambiguous=ambiguous,
        opened=opened, candidates=candidates, old_bounds=old_bounds,
        new_bounds=new_bounds,
    )
    if change_log is not None:
        change_log.metadata["target"] = dict(info)

    model.c[id_biomass] = 0
    model.c[id_target] = 1

    success, x, _ = _maximize(model.c, model.S, model.b,
                              np.zeros((0, model.S.shape[1])),
                              np.zeros(0), model.lb, model.ub,
                              backend=backend)
    if not success:
        model.c[id_biomass] = 1
        model.c[id_target] = 0
        return (model, id_target, 0.0, info) if return_info else (model, id_target, 0.0)

    TMPR = float(x[id_target])
    model.c[id_biomass] = 1
    model.c[id_target] = 0
    return (model, id_target, TMPR, info) if return_info else (model, id_target, TMPR)


def _unique_reaction_id(model: Model, base_id: str) -> str:
    """Return ``base_id`` or a suffixed variant that is not already used."""
    existing = {str(name) for name in model.rxns}
    if base_id not in existing:
        return base_id
    suffix = 2
    while f"{base_id}_{suffix}" in existing:
        suffix += 1
    return f"{base_id}_{suffix}"


def _add_target_reaction(model: Model, id_met: int, kind: str) -> int:
    """Append a virtual target outlet reaction and return its index."""
    secretion = kind == "secretion"
    id_target = len(model.rxns)
    base_id = "EX_target" if secretion else "DM_target"
    reaction_id = _unique_reaction_id(model, base_id)
    model.S = np.hstack((model.S, np.zeros((model.S.shape[0], 1))))
    # In both cases a positive flux represents production of the target: an
    # exchange/demand reaction consumes the metabolite from the network.
    model.S[id_met, id_target] = -1
    model.rxns = list(model.rxns) + [reaction_id]
    model.grRules = list(model.grRules) + [""]
    if len(model.rxnNames) > 0:
        model.rxnNames = list(model.rxnNames) + [
            "Virtual target secretion outlet" if secretion
            else "Virtual target intracellular demand"
        ]
    if len(model.subSystems) > 0:
        model.subSystems = list(model.subSystems) + [""]
    if secretion:
        model.lb = np.append(model.lb, -1000.0)
        model.ub = np.append(model.ub, 1000.0)
    else:
        # A demand reaction is an irreversible intracellular sink.
        model.lb = np.append(model.lb, 0.0)
        model.ub = np.append(model.ub, 1000.0)
    model.c = np.append(model.c, 0.0)
    if model.rev is not None:
        model.rev = np.append(np.asarray(model.rev).ravel(), 1.0 if secretion else 0.0)
    if model.rxnGeneMat is not None:
        model.rxnGeneMat = np.vstack(
            (model.rxnGeneMat, np.zeros((1, model.rxnGeneMat.shape[1])))
        )
    return id_target
