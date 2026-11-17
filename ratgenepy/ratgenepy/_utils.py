"""Small helpers for MATLAB data conversion, matrix building and indexing."""

from __future__ import annotations

from typing import List, Sequence

import numpy as np
from scipy import sparse

def _as_str(value) -> str:
    """Convert a MATLAB char/cell/string loaded by scipy into a Python str."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ""
        return _as_str(value.ravel()[0])
    if hasattr(value, "ravel"):
        try:
            arr = np.asarray(value).ravel()
            if arr.size == 0:
                return ""
            return _as_str(arr[0])
        except Exception:
            pass
    return str(value)


def _as_cellstr(value) -> List[str]:
    """Convert a MATLAB cell array of strings into a list of Python strings."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return []
        return [_as_str(item) for item in value.ravel()]
    if isinstance(value, (list, tuple)):
        return [_as_str(item) for item in value]
    return [_as_str(value)]


def _as_vec(value, dtype=float) -> np.ndarray:
    """Convert a MATLAB numeric array into a flat numpy vector."""
    if value is None:
        return np.array([], dtype=dtype)
    if sparse.issparse(value):
        value = value.toarray()
    return np.asarray(value, dtype=dtype).ravel()


def _as_matrix(value, dtype=float) -> np.ndarray:
    """Convert a MATLAB numeric matrix (possibly sparse) into a dense matrix."""
    if value is None:
        return np.zeros((0, 0), dtype=dtype)
    if sparse.issparse(value):
        value = value.toarray()
    matrix = np.asarray(value, dtype=dtype)
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    return matrix


def _normalize_vtype(vtype, n: int) -> List[str]:
    if vtype is None:
        return ["C"] * n
    if isinstance(vtype, str):
        if len(vtype) == n:
            return list(vtype)
        if len(vtype) == 1:
            return [vtype] * n
        raise ValueError("vtype length does not match the number of variables")
    vtype = list(vtype)
    if len(vtype) != n:
        raise ValueError("vtype length does not match the number of variables")
    return vtype


def _blkdiag(matrix: np.ndarray, block: np.ndarray) -> np.ndarray:
    """Equivalent of MATLAB blkdiag for two dense matrices."""
    matrix = np.asarray(matrix, dtype=float)
    block = np.asarray(block, dtype=float)
    if matrix.shape[0] == 0 and matrix.shape[1] == 0:
        return block.copy()
    if block.shape[0] == 0 and block.shape[1] == 0:
        return matrix.copy()
    out = np.zeros((matrix.shape[0] + block.shape[0],
                    matrix.shape[1] + block.shape[1]), dtype=float)
    out[:matrix.shape[0], :matrix.shape[1]] = matrix
    out[matrix.shape[0]:, matrix.shape[1]:] = block
    return out


def _vstack(matrix: np.ndarray, block: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    block = np.asarray(block, dtype=float)
    if matrix.shape[0] == 0 and matrix.shape[1] == 0:
        return block.copy()
    if block.shape[0] == 0 and block.shape[1] == 0:
        return matrix.copy()
    if matrix.shape[1] != block.shape[1]:
        raise ValueError("inconsistent number of columns")
    return np.vstack((matrix, block))


def _gene_index(genes: Sequence[str], name: str) -> int:
    for i, gene in enumerate(genes):
        if gene == name:
            return i
    for i, gene in enumerate(genes):
        if name in gene:
            return i
    raise ValueError(f"gene {name!r} is not in the model")


def _first_match(values: Sequence[str], name: str) -> int:
    matches = [i for i, value in enumerate(values) if value == name]
    if not matches:
        raise NameError(f"{name!r} is not in the model")
    return matches[0]
