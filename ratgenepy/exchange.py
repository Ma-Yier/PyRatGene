"""Detection or introduction of the target metabolite exchange reaction."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .model import Model
from .solvers import _maximize

def introExchange(input_model: Model, id_biomass: int, id_input: Sequence[int],
                  id_met: int, backend: str = "scip"):
    """Detect or add the exchange reaction of the target metabolite.

    Returns ``(model, id_target, TMPR)``.
    """
    model = input_model.copy()
    id_input = list(np.asarray(id_input).ravel())

    if len(id_input) == 1:
        id_carbon = int(id_input[0])
        id_oxygen = -1
    elif len(id_input) == 2:
        id_carbon = int(id_input[0])
        id_oxygen = int(id_input[1])
    else:
        raise ValueError("Input id error. At most 2 allowed.")

    name_met = model.mets[id_met]
    ex_rxn = "EX_" + str(name_met)
    dx_rxn = "DM_" + str(name_met)

    id_target = None
    if ex_rxn in model.rxns:
        id_target = model.rxns.index(ex_rxn)
    if dx_rxn in model.rxns:
        id_target = model.rxns.index(dx_rxn)

    if id_target is not None:
        if id_target == id_carbon or id_target == id_oxygen:
            raise ValueError("target metabolite cannot be the carbon or oxygen source")
    else:
        id_target = len(model.rxns)
        model.S = np.hstack((model.S, np.zeros((model.S.shape[0], 1))))
        model.S[id_met, id_target] = -1
        model.rxns = list(model.rxns) + ["EX_target"]
        model.grRules = list(model.grRules) + [""]
        if len(model.rxnNames) > 0:
            model.rxnNames = list(model.rxnNames) + ["EX_target"]
        if len(model.subSystems) > 0:
            model.subSystems = list(model.subSystems) + [""]
        model.lb = np.append(model.lb, -1000.0)
        model.ub = np.append(model.ub, 1000.0)
        model.c = np.append(model.c, 0.0)
        if model.rev is not None:
            model.rev = np.append(np.asarray(model.rev).ravel(), 1.0)
        if model.rxnGeneMat is not None:
            model.rxnGeneMat = np.vstack(
                (model.rxnGeneMat, np.zeros((1, model.rxnGeneMat.shape[1])))
            )

    model.c[id_biomass] = 0
    model.c[id_target] = 1

    success, x, _ = _maximize(model.c, model.S, model.b,
                              np.zeros((0, model.S.shape[1])),
                              np.zeros(0), model.lb, model.ub,
                              backend=backend)
    if not success:
        model.c[id_biomass] = 1
        model.c[id_target] = 0
        return model, id_target, 0.0

    TMPR = float(x[id_target])
    model.c[id_biomass] = 1
    model.c[id_target] = 0
    return model, id_target, TMPR
