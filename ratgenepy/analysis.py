"""Worst-case and best-case target-rate analysis for a strategy."""

from __future__ import annotations

import numpy as np

from ._utils import _as_vec
from .gpr import verifyRatioGene
from .model import Model
from .solvers import _maximize, _minimize

def boundAnalysis(model: Model, kogene, id_biomass: int, id_target: int,
                  backend: str = "scip"):
    """Return the worst-case and best-case target rates for a strategy."""
    minTarget = 0.0
    maxTarget = 0.0

    kogene = _as_vec(kogene)
    var_x = verifyRatioGene(model, kogene)
    model1 = model.copy()

    index = np.where(var_x == 1)[0]
    model1.lb[index] = 0
    model1.ub[index] = 0
    model1.lb[id_biomass] = 0

    success, _, biomass = _maximize(
        model1.c, model1.S, model1.b, np.zeros((0, model1.S.shape[1])),
        np.zeros(0), model1.lb, model1.ub, backend=backend
    )
    if not success:
        return minTarget, maxTarget

    model2 = model1.copy()
    model2.lb[id_biomass] = biomass
    model2.ub[id_biomass] = biomass
    model2.c[id_biomass] = 0
    model2.c[id_target] = 1

    _, _, maxTarget = _maximize(
        model2.c, model2.S, model2.b, np.zeros((0, model2.S.shape[1])),
        np.zeros(0), model2.lb, model2.ub, backend=backend
    )
    if maxTarget <= 0.001:
        maxTarget = 0.0

    _, _, minTarget = _minimize(
        model2.c, model2.S, model2.b, np.zeros((0, model2.S.shape[1])),
        np.zeros(0), model2.lb, model2.ub, backend=backend
    )
    if minTarget <= 0.001:
        minTarget = 0.0

    return minTarget, maxTarget
