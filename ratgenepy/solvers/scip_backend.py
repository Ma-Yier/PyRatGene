"""SCIP backend for LP/MILP problems."""

from __future__ import annotations

import weakref
from typing import Any, List, Tuple

import numpy as np

from ._common import normalize_problem

_SCIP_VARIABLES = weakref.WeakKeyDictionary()


def _require_scip():
    try:
        from pyscipopt import Model as ScipModel
        from pyscipopt import quicksum
    except ImportError as exc:
        raise ImportError(
            "PySCIPOpt is not installed. Install it with `pip install pyscipopt` "
            "or choose another backend with backend='cplex' or backend='gurobi'."
        ) from exc
    return ScipModel, quicksum


def _make_problem(objfunc, Aeq, beq, A, b, lb, ub, vtype,
                  timeLimit=None, pool=None) -> Tuple[Any, list]:
    ScipModel, quicksum = _require_scip()
    n = len(lb)
    scip = ScipModel()
    scip.hideOutput()

    variables = []
    for i in range(n):
        variables.append(
            scip.addVar(name=f"x_{i}", vtype=vtype[i], lb=float(lb[i]),
                        ub=float(ub[i]))
        )
    # SCIP may reorder its variables internally; keep the construction order so
    # that solution vectors always follow [flux; indicators; ...].
    _SCIP_VARIABLES[scip] = variables

    Aeq = np.asarray(Aeq, dtype=float)
    beq = np.asarray(beq, dtype=float).ravel()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).ravel()

    for i in range(Aeq.shape[0]):
        if np.all(Aeq[i, :] == 0):
            continue
        expr = quicksum(float(Aeq[i, j]) * variables[j]
                        for j in range(n) if Aeq[i, j] != 0)
        scip.addCons(expr == float(beq[i]), name=f"eq_{i}")

    for i in range(A.shape[0]):
        if np.all(A[i, :] == 0):
            continue
        expr = quicksum(float(A[i, j]) * variables[j]
                        for j in range(n) if A[i, j] != 0)
        scip.addCons(expr <= float(b[i]), name=f"ineq_{i}")

    objfunc = np.asarray(objfunc, dtype=float).ravel()
    scip.setObjective(quicksum(float(objfunc[j]) * variables[j]
                               for j in range(n) if objfunc[j] != 0),
                      sense="minimize")

    if timeLimit is not None and np.isfinite(timeLimit):
        scip.setParam("limits/time", float(timeLimit))

    return scip, variables


def _solution_vector(scip, variables, sol) -> np.ndarray:
    return np.array([scip.getSolVal(sol, var) for var in variables], dtype=float)


def scip_solver(objfunc, **kwargs):
    """Solve an LP or MILP with SCIP.

    Returns ``(scip_model, x)``.  See :func:`ratgenepy.solvers.solver` for the
    accepted keyword arguments.
    """
    problem = normalize_problem(objfunc, kwargs)
    scip, variables = _make_problem(
        problem["obj"],
        problem["Aeq"], problem["beq"],
        problem["A"], problem["b"],
        problem["lb"], problem["ub"],
        problem["vtype"],
        timeLimit=problem["timeLimit"],
        pool=problem["pool"],
    )
    scip.optimize()

    if scip.getNSols() == 0:
        return scip, None

    best = scip.getBestSol()
    return scip, _solution_vector(scip, variables, best)


def pool_solutions(scip, variables=None, limit: int = 10) -> List[np.ndarray]:
    """Return up to ``limit`` stored SCIP solutions, best objective first."""
    if scip.getNSols() == 0:
        return []

    if variables is None:
        variables = _SCIP_VARIABLES.get(scip, scip.getVars())

    def objective(sol):
        try:
            return scip.getSolObjVal(sol)
        except Exception:
            return float("inf")

    solutions = sorted(scip.getSols(), key=objective)
    return [_solution_vector(scip, variables, sol) for sol in solutions[:limit]]


def _make_scip_variables(scip):
    """Return the SCIP variables in construction order."""
    return _SCIP_VARIABLES.get(scip, scip.getVars())
