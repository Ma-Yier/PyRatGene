"""Gurobi backend for LP/MILP problems.

The Gurobi Python API is optional.  Importing this module does not import
``gurobipy``; an informative error is raised only when the backend is used and
the API is unavailable.
"""

from __future__ import annotations

from typing import List

import numpy as np

from ._common import normalize_problem


def _require_gurobi():
    try:
        import gurobipy as gp
        from gurobipy import GRB
    except ImportError as exc:
        raise ImportError(
            "Gurobi is not installed. Install `gurobipy` and a valid Gurobi "
            "license, or choose another backend with backend='scip' or "
            "backend='cplex'."
        ) from exc
    return gp, GRB


def _bound(GRB, value: float) -> float:
    if np.isposinf(value):
        return float(GRB.INFINITY)
    if np.isneginf(value):
        return float(-GRB.INFINITY)
    return float(value)


def _add_constraints(model, gp, variables, matrix: np.ndarray,
                     rhs: np.ndarray, sense: str) -> None:
    for i in range(matrix.shape[0]):
        row = matrix[i]
        indices = np.nonzero(row)[0]
        if indices.size == 0:
            continue
        expr = gp.quicksum(float(row[j]) * variables[int(j)]
                           for j in indices)
        if sense == "E":
            model.addConstr(expr == float(rhs[i]), name=f"eq_{i}")
        else:
            model.addConstr(expr <= float(rhs[i]), name=f"ineq_{i}")


def gurobi_solver(objfunc, **kwargs):
    """Solve an LP or MILP with Gurobi.

    Returns ``(gurobi_model, x)``.  See
    :func:`ratgenepy.solvers.solver` for the accepted keyword arguments.
    """
    gp, GRB = _require_gurobi()
    problem = normalize_problem(objfunc, kwargs)

    model = gp.Model()
    model.Params.OutputFlag = 0

    type_map = {
        "C": GRB.CONTINUOUS,
        "B": GRB.BINARY,
        "I": GRB.INTEGER,
    }
    vtypes = [
        type_map.get(vtype, GRB.CONTINUOUS)
        for vtype in problem["vtype"]
    ]

    variables = model.addVars(
        list(range(len(problem["lb"]))),
        lb=[_bound(GRB, value) for value in problem["lb"]],
        ub=[_bound(GRB, value) for value in problem["ub"]],
        vtype=vtypes,
        name="x",
    )
    variable_list = [variables[i] for i in range(len(problem["lb"]))]

    model.setObjective(
        gp.quicksum(
            float(problem["obj"][j]) * variable_list[j]
            for j in range(len(variable_list))
            if problem["obj"][j] != 0
        ),
        GRB.MINIMIZE,
    )

    _add_constraints(
        model, gp, variable_list, problem["Aeq"], problem["beq"], "E"
    )
    _add_constraints(
        model, gp, variable_list, problem["A"], problem["b"], "L"
    )

    if (problem["timeLimit"] is not None
            and np.isfinite(problem["timeLimit"])):
        model.Params.TimeLimit = float(problem["timeLimit"])

    if problem["opt"] == "milp" and problem["pool"] > 1:
        model.Params.PoolSearchMode = 1
        model.Params.PoolSolutions = int(problem["pool"])

    model.optimize()

    if model.SolCount == 0:
        # Some continuous models report SolCount == 0 even after a successful
        # solve; fall back to the current solution when one is available.
        try:
            return model, np.asarray(model.X, dtype=float)
        except Exception:
            return model, None

    return model, np.asarray(model.X, dtype=float)


def pool_solutions(model, limit: int = 10) -> List[np.ndarray]:
    """Return up to ``limit`` stored Gurobi solutions, best objective first."""
    if model.SolCount == 0:
        try:
            return [np.asarray(model.X, dtype=float)][:limit]
        except Exception:
            return []

    variables = model.getVars()
    if model.SolCount == 1:
        return [
            np.asarray([variable.X for variable in variables], dtype=float)
        ][:limit]

    values = []
    objectives = []
    original = model.Params.SolutionNumber
    try:
        for i in range(model.SolCount):
            model.Params.SolutionNumber = i
            values.append(
                np.asarray([variable.Xn for variable in variables], dtype=float)
            )
            try:
                # Gurobi 13 renamed PoolObjVal to PoolNObjVal; keep both paths.
                try:
                    objective = model.PoolNObjVal
                except Exception:
                    objective = model.PoolObjVal
                objectives.append(float(objective))
            except Exception:
                objectives.append(float("inf"))
    finally:
        model.Params.SolutionNumber = original

    order = np.argsort(objectives)
    return [values[int(i)] for i in order[:limit]]
