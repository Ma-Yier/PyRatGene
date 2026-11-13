"""Multi-backend LP/MILP solver interface.

Supported backends:

- ``scip``   : PySCIPOpt (default)
- ``cplex``  : IBM ILOG CPLEX Python API
- ``gurobi`` : Gurobi Python API (``gurobipy``)

Use :func:`solver` with ``backend="scip"``, ``"cplex"`` or ``"gurobi"``.
The backend-specific functions :func:`scip_solver`, :func:`cplex_solver` and
:func:`gurobi_solver` are also available directly.
"""

from __future__ import annotations

import importlib.util
from typing import List

import numpy as np

from ..model import Model
from .cplex_backend import cplex_solver
from .cplex_backend import pool_solutions as _cplex_pool_solutions
from .gurobi_backend import gurobi_solver
from .gurobi_backend import pool_solutions as _gurobi_pool_solutions
from .scip_backend import scip_solver
from .scip_backend import pool_solutions as _scip_pool_solutions

_BACKENDS = {
    "scip": scip_solver,
    "cplex": cplex_solver,
    "gurobi": gurobi_solver,
}

_BACKEND_MODULES = {
    "scip": "pyscipopt",
    "cplex": "cplex",
    "gurobi": "gurobipy",
}


def available_solvers() -> List[str]:
    """Return the solver backends that are currently importable."""
    available = []
    for name, module_name in _BACKEND_MODULES.items():
        try:
            if importlib.util.find_spec(module_name) is not None:
                available.append(name)
        except (ImportError, ValueError):
            continue
    return available


def solver(objfunc, backend="scip", **kwargs):
    """Solve an LP or MILP problem.

    Parameters
    ----------
    objfunc : array-like
        Linear objective vector.  The objective is always minimized.
    backend : {"scip", "cplex", "gurobi"}
        Solver backend.  Default is ``"scip"``.
    **kwargs
        ``Aeq, beq`` : equality constraints ``Aeq x == beq``
        ``A, b``     : inequality constraints ``A x <= b``
        ``lb, ub``   : variable bounds
        ``opt``      : ``"lp"`` or ``"milp"``
        ``vtype``    : variable types for MILP (``"C"``, ``"B"`` or ``"I"``)
        ``pool``     : maximum number of candidate solutions kept for verification
        ``timeLimit``: solver time limit in seconds
        ``integerTolerance``: integrality tolerance for MILP.  If omitted,
            the smallest value supported by the selected backend is used:
            SCIP ``numerics/feastol`` = ``1e-17``, CPLEX
            ``mip.tolerances.integrality`` = ``0.0``, Gurobi
            ``IntFeasTol`` = ``1e-9``.

    Returns
    -------
    (handle, x)
        ``handle`` is the backend model object and ``x`` is the best solution
        vector, or ``None`` when the problem has no feasible solution.
    """
    key = str(backend).lower()
    if key not in _BACKENDS:
        raise ValueError(
            f"unknown solver backend {backend!r}; choose from "
            f"{', '.join(sorted(_BACKENDS))}"
        )
    return _BACKENDS[key](objfunc, **kwargs)


def pool_solutions(handle, limit: int = 10) -> List[np.ndarray]:
    """Return up to ``limit`` stored solutions, best objective first."""
    module = type(handle).__module__.lower()

    if module.startswith("cplex"):
        return _cplex_pool_solutions(handle, limit=limit)
    if module.startswith("gurobipy"):
        return _gurobi_pool_solutions(handle, limit=limit)
    if module.startswith("pyscipopt"):
        return _scip_pool_solutions(handle, limit=limit)

    # Fallback for backend objects exposed through a wrapper module.
    if hasattr(handle, "SolCount"):
        return _gurobi_pool_solutions(handle, limit=limit)
    if hasattr(handle, "solution") and hasattr(handle.solution, "pool"):
        return _cplex_pool_solutions(handle, limit=limit)
    if hasattr(handle, "getNSols") and hasattr(handle, "getBestSol"):
        return _scip_pool_solutions(handle, limit=limit)

    raise TypeError(f"unsupported solver handle: {type(handle)!r}")


def _maximize(c, Aeq, beq, A, b, lb, ub, timeLimit=None, backend="scip"):
    """Maximize c*x; return ``(success, x, objective)``."""
    _, x = solver(
        -np.asarray(c, dtype=float).ravel(),
        backend=backend,
        Aeq=Aeq, beq=beq, A=A, b=b, lb=lb, ub=ub,
        opt="lp", timeLimit=timeLimit,
    )
    if x is None:
        return False, None, 0.0
    return True, x, float(np.dot(np.asarray(c).ravel(), x))


def _minimize(c, Aeq, beq, A, b, lb, ub, timeLimit=None, backend="scip"):
    """Minimize c*x; return ``(success, x, objective)``."""
    _, x = solver(
        np.asarray(c, dtype=float).ravel(),
        backend=backend,
        Aeq=Aeq, beq=beq, A=A, b=b, lb=lb, ub=ub,
        opt="lp", timeLimit=timeLimit,
    )
    if x is None:
        return False, None, 0.0
    return True, x, float(np.dot(np.asarray(c).ravel(), x))


def optimizeCbModel(model: Model, backend="scip"):
    """Maximize the model objective; return a small result object."""

    class Result:
        pass

    result = Result()
    success, x, value = _maximize(
        model.c, model.S, model.b,
        np.zeros((0, model.S.shape[1])),
        np.zeros(0), model.lb, model.ub,
        backend=backend,
    )
    result.stat = 1 if success else 0
    result.x = x if x is not None else np.zeros(len(model.rxns))
    result.f = value if success else 0.0
    return result


__all__ = [
    "available_solvers",
    "cplex_solver",
    "gurobi_solver",
    "optimizeCbModel",
    "pool_solutions",
    "scip_solver",
    "solver",
]
