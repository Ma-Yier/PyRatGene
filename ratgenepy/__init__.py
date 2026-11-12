"""RatGene: growth-to-production ratio-based design algorithms.

The package mirrors the MATLAB RatGene implementation.  The main entry point
is :func:`ratgenepy.core.RatGene`; the remaining modules expose the individual
steps of the algorithm (GPR parsing, ratio-based search, strategy reduction,
model integration and bounds analysis).
"""

from __future__ import annotations

from .analysis import boundAnalysis
from .candidates import getCandidateAdd, getCandidateDel
from .core import RatGene
from .exchange import introExchange
from .gpr import constructMatrix, locateGene, maxParseGPR, verifyRatioGene
from .integrate import findDiffGene, uniteGeneModel
from .model import Model
from .ratio import RatMethodRxn, ratioGene, ratioMethod, verifyRxn
from .reduction import rmRedundancy
from .solvers import (
    available_solvers,
    cplex_solver,
    gurobi_solver,
    optimizeCbModel,
    pool_solutions,
    scip_solver,
    solver,
)
from .verification import verifyGeneKnock, verifyKnock

__version__ = "1.0.0"

__all__ = [
    "Model",
    "RatGene",
    "RatMethodRxn",
    "available_solvers",
    "boundAnalysis",
    "constructMatrix",
    "cplex_solver",
    "findDiffGene",
    "getCandidateAdd",
    "getCandidateDel",
    "gurobi_solver",
    "introExchange",
    "locateGene",
    "maxParseGPR",
    "optimizeCbModel",
    "pool_solutions",
    "ratioGene",
    "ratioMethod",
    "rmRedundancy",
    "scip_solver",
    "solver",
    "uniteGeneModel",
    "verifyGeneKnock",
    "verifyKnock",
    "verifyRatioGene",
    "verifyRxn",
]
