"""Shared argument parsing for LP/MILP solver backends."""

from __future__ import annotations

from typing import Any, Dict

import numpy as np

from .._utils import _as_vec, _normalize_vtype


# Smallest/largest values accepted by each backend for the integrality
# (or closest equivalent) tolerance.  CPLEX documents [0.0, 0.5] for
# mip.tolerances.integrality, Gurobi documents [1e-9, 1e-1] for IntFeasTol,
# and SCIP uses numerics/feastol with range [1e-17, 1e-3].
_MIN_INTEGER_TOLERANCE = {
    "scip": 1e-10,
    "cplex": 1e-12,
    "gurobi": 1e-9,
}
_MAX_INTEGER_TOLERANCE = {
    "scip": 1e-3,
    "cplex": 0.5,
    "gurobi": 1e-1,
}


def resolve_integer_tolerance(problem: Dict[str, Any], backend: str) -> float:
    """Return the integrality tolerance, clamped to the backend's range.

    ``integerTolerance=None`` means "use the smallest value supported by the
    selected backend".
    """
    lower = _MIN_INTEGER_TOLERANCE[backend]
    upper = _MAX_INTEGER_TOLERANCE[backend]
    value = problem.get("integerTolerance")
    if value is None:
        return lower
    return float(min(max(float(value), lower), upper))


def _as_2d(value) -> np.ndarray:
    """Convert a constraint matrix to a dense 2-D array."""
    if value is None:
        return np.zeros((0, 0), dtype=float)
    matrix = np.asarray(value, dtype=float)
    if matrix.size == 0:
        return np.zeros((0, 0), dtype=float)
    if matrix.ndim == 0:
        matrix = matrix.reshape(1, 1)
    elif matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    return matrix


def normalize_problem(objfunc, kwargs) -> Dict[str, Any]:
    """Normalize the common solver arguments used by all backends."""
    opt = kwargs.get("opt")
    if opt is None:
        raise KeyError("should assign optimization type LP or MILP")
    opt = str(opt).lower()
    if opt not in ("lp", "milp"):
        raise ValueError("opt should be either lp or milp")

    Aeq = _as_2d(kwargs.get("Aeq", np.zeros((0, 0))))
    beq = _as_vec(kwargs.get("beq", np.zeros(0)))
    A = _as_2d(kwargs.get("A", np.zeros((0, 0))))
    b = _as_vec(kwargs.get("b", np.zeros(0)))
    lb = _as_vec(kwargs.get("lb"))
    ub = _as_vec(kwargs.get("ub"))

    n = len(lb)
    if len(ub) != n:
        raise ValueError("lb and ub must have the same length")

    objective = np.asarray(objfunc, dtype=float).ravel()
    if objective.size != n:
        raise ValueError("objective vector length does not match lb/ub")

    if Aeq.shape[0] != beq.size:
        raise ValueError("Aeq and beq have inconsistent row counts")
    if A.shape[0] != b.size:
        raise ValueError("A and b have inconsistent row counts")
    if Aeq.shape[0] and Aeq.shape[1] != n:
        raise ValueError("Aeq has the wrong number of columns")
    if A.shape[0] and A.shape[1] != n:
        raise ValueError("A has the wrong number of columns")

    if opt == "lp":
        vtype = ["C"] * n
    else:
        vtype = [str(v).upper() for v in _normalize_vtype(kwargs.get("vtype"), n)]

    integer_tolerance = kwargs.get("integerTolerance")
    if integer_tolerance is not None:
        integer_tolerance = float(integer_tolerance)
        if not np.isfinite(integer_tolerance) or integer_tolerance < 0:
            raise ValueError("integerTolerance must be a finite non-negative value")

    return {
        "obj": objective,
        "Aeq": Aeq,
        "beq": beq,
        "A": A,
        "b": b,
        "lb": lb,
        "ub": ub,
        "vtype": vtype,
        "opt": opt,
        "timeLimit": kwargs.get("timeLimit"),
        "pool": max(0, int(kwargs.get("pool", 0) or 0)),
        "integerTolerance": integer_tolerance,
    }
