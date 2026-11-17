"""CPLEX backend for LP/MILP problems.

The CPLEX Python API is optional.  Importing this module does not import
``cplex``; an informative error is raised only when the backend is used and the
API is unavailable.
"""

from __future__ import annotations

from typing import List

import numpy as np

from ._common import normalize_problem, resolve_integer_tolerance


def _require_cplex():
    try:
        import cplex
    except ImportError as exc:
        raise ImportError(
            "The CPLEX Python API is not installed. Install IBM ILOG CPLEX "
            "Optimization Studio and make its Python API importable as "
            "`cplex`, or choose another backend with backend='scip' or "
            "backend='gurobi'."
        ) from exc
    return cplex


def _bound(cplex, value: float) -> float:
    if np.isposinf(value):
        return float(cplex.infinity)
    if np.isneginf(value):
        return float(-cplex.infinity)
    return float(value)


def _add_constraints(cpx, cplex, matrix: np.ndarray, rhs: np.ndarray,
                     sense: str) -> None:
    lin_expr = []
    senses = []
    rhs_values = []

    for i in range(matrix.shape[0]):
        row = matrix[i]
        indices = np.nonzero(row)[0]
        if indices.size == 0:
            continue
        lin_expr.append(
            cplex.SparsePair(
                ind=[int(j) for j in indices],
                val=[float(row[j]) for j in indices],
            )
        )
        senses.append(sense)
        rhs_values.append(float(rhs[i]))

    if lin_expr:
        cpx.linear_constraints.add(
            lin_expr=lin_expr,
            senses="".join(senses),
            rhs=rhs_values,
        )


def _current_solution(cpx):
    try:
        if cpx.solution.is_primal_feasible():
            return np.asarray(cpx.solution.get_values(), dtype=float)
    except Exception:
        pass
    try:
        return np.asarray(cpx.solution.get_values(), dtype=float)
    except Exception:
        return None


def _pool_values(cpx) -> List[np.ndarray]:
    try:
        count = cpx.solution.pool.get_num()
    except Exception:
        return []

    if count <= 0:
        return []

    values = []
    objectives = []
    for i in range(count):
        try:
            values.append(
                np.asarray(cpx.solution.pool.get_values(i), dtype=float)
            )
        except Exception:
            continue
        try:
            objectives.append(float(cpx.solution.pool.get_objective_value(i)))
        except Exception:
            objectives.append(float("inf"))

    if not values:
        return []
    if len(objectives) != len(values):
        objectives = [float("inf")] * len(values)

    order = np.argsort(objectives)
    return [values[int(i)] for i in order]


def cplex_solver(objfunc, **kwargs):
    """Solve an LP or MILP with CPLEX.

    Returns ``(cplex_model, x)``.  See
    :func:`ratgenepy.solvers.solver` for the accepted keyword arguments.
    """
    cplex = _require_cplex()
    problem = normalize_problem(objfunc, kwargs)

    cpx = cplex.Cplex()
    try:
        cpx.set_log_stream(None)
        cpx.set_error_stream(None)
        cpx.set_warning_stream(None)
    except Exception:
        pass

    cpx.objective.set_sense(cpx.objective.sense.minimize)

    if problem["opt"] == "milp":
        cpx.parameters.mip.tolerances.integrality.set(
            resolve_integer_tolerance(problem, "cplex")
        )

    type_map = {
        "C": cpx.variables.type.continuous,
        "B": cpx.variables.type.binary,
        "I": cpx.variables.type.integer,
    }
    types = [
        type_map.get(vtype, cpx.variables.type.continuous)
        for vtype in problem["vtype"]
    ]

    cpx.variables.add(
        obj=[float(value) for value in problem["obj"]],
        lb=[_bound(cplex, value) for value in problem["lb"]],
        ub=[_bound(cplex, value) for value in problem["ub"]],
        types=types,
    )

    _add_constraints(cpx, cplex, problem["Aeq"], problem["beq"], "E")
    _add_constraints(cpx, cplex, problem["A"], problem["b"], "L")

    if (problem["timeLimit"] is not None
            and np.isfinite(problem["timeLimit"])):
        cpx.parameters.timelimit.set(float(problem["timeLimit"]))

    pool_size = problem["pool"]
    if problem["opt"] == "milp" and pool_size > 1:
        try:
            cpx.parameters.mip.pool.capacity.set(int(pool_size))
        except Exception:
            pass

    use_pool = (
        problem["opt"] == "milp"
        and pool_size > 1
        and hasattr(cpx, "populate_solution_pool")
    )
    if use_pool:
        cpx.populate_solution_pool()
    else:
        cpx.solve()

    values = _pool_values(cpx)
    if values:
        return cpx, values[0]

    return cpx, _current_solution(cpx)


def pool_solutions(cpx, limit: int = 10) -> List[np.ndarray]:
    """Return up to ``limit`` stored CPLEX solutions, best objective first."""
    values = _pool_values(cpx)
    if not values:
        current = _current_solution(cpx)
        if current is not None:
            values = [current]
    return values[:limit]
