"""Verification of gene and reaction modification strategies."""

from __future__ import annotations

import numpy as np

from ._utils import _as_vec
from .gpr import verifyRatioGene
from .model import Model
from .solvers import _maximize

def verifyGeneKnock(model: Model, input_x, indGPR, id_biomass: int,
                    id_target: int, min_bound: float,
                    backend: str = "scip"):
    """Validate a gene deletion/addition strategy.

    Returns ``(x_target, gene_strategy)``.
    """
    work = model.copy()
    input_x = _as_vec(input_x)
    input_xg = np.zeros_like(input_x)
    input_xg[input_x > 0.1] = 1

    var_x = verifyRatioGene(work, input_xg)
    index = np.where(var_x == 1)[0]
    work.lb[index] = 0
    work.ub[index] = 0
    work.lb[id_biomass] = 0

    success, x, value = _maximize(work.c, work.S, work.b,
                                  np.zeros((0, work.S.shape[1])),
                                  np.zeros(0), work.lb, work.ub,
                                  backend=backend)
    if success and value >= min_bound:
        return float(x[id_target]), input_xg
    return 0.0, input_xg


def verifyKnock(model: Model, input_x, id_biomass: int, id_target: int,
                min_bound: float, backend: str = "scip"):
    """Validate a reaction deletion/addition strategy.

    Returns ``(x_target, knockout)``.
    """
    work = model.copy()
    var_x = _as_vec(input_x).copy()
    var_x[np.abs(var_x) < 1e-7] = 0

    work.lb[var_x == 0] = 0
    work.ub[var_x == 0] = 0
    work.lb[id_biomass] = 0

    knockout = np.ones_like(var_x)
    knockout[var_x == 0] = 0

    success, x, value = _maximize(work.c, work.S, work.b,
                                  np.zeros((0, work.S.shape[1])),
                                  np.zeros(0), work.lb, work.ub,
                                  backend=backend)
    if success and value > min_bound:
        return float(x[id_target]), knockout
    return 0.0, knockout
