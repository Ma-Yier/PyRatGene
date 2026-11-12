"""Backward-compatible SCIP interface.

New code should use :mod:`ratgenepy.solvers`, which supports SCIP, CPLEX and
Gurobi through a common ``solver`` function.  This module keeps the original
SCIP-specific names available for older code.
"""

from __future__ import annotations

from .solvers import (
    _maximize,
    _minimize,
    available_solvers,
    cplex_solver,
    gurobi_solver,
    optimizeCbModel,
    pool_solutions,
    solver,
)
from .solvers.scip_backend import (
    _SCIP_VARIABLES,
    _make_problem,
    _make_scip_variables,
    _solution_vector,
    scip_solver,
)


def _pool_solutions(scip, variables=None, limit: int = 10):
    """Compatibility wrapper for the old SCIP solution-pool helper."""
    return pool_solutions(scip, limit=limit)


__all__ = [
    "_SCIP_VARIABLES",
    "_make_problem",
    "_make_scip_variables",
    "_maximize",
    "_minimize",
    "_pool_solutions",
    "_solution_vector",
    "available_solvers",
    "cplex_solver",
    "gurobi_solver",
    "optimizeCbModel",
    "pool_solutions",
    "scip_solver",
    "solver",
]
