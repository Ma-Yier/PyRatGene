"""COBRA-style model container and .mat model loading."""

from __future__ import annotations

import copy
from typing import Optional

import numpy as np
import scipy.io as sco

try:  # scipy >= 1.8
    from scipy.io.matlab import mat_struct
except Exception:  # pragma: no cover - very old scipy
    mat_struct = None

from ._utils import _as_cellstr, _as_matrix, _as_str, _as_vec

class Model:
    """A small COBRA-style model container used by the RatGene functions."""

    def __init__(self, **kwargs):
        self.S = kwargs.get("S", np.zeros((0, 0)))
        self.b = kwargs.get("b", np.zeros(self.S.shape[0]))
        self.c = kwargs.get("c", np.zeros(self.S.shape[1]))
        self.lb = kwargs.get("lb", np.zeros(self.S.shape[1]))
        self.ub = kwargs.get("ub", np.zeros(self.S.shape[1]))
        self.mets = kwargs.get("mets", [])
        self.metNames = kwargs.get("metNames", [])
        self.metFormulas = kwargs.get("metFormulas", [])
        self.metCharge = kwargs.get("metCharge", None)
        self.rxns = kwargs.get("rxns", [])
        self.rxnNames = kwargs.get("rxnNames", [])
        self.subSystems = kwargs.get("subSystems", [])
        self.grRules = kwargs.get("grRules", [])
        self.genes = kwargs.get("genes", [])
        self.rev = kwargs.get("rev", None)
        self.rxnGeneMat = kwargs.get("rxnGeneMat", None)
        self.description = kwargs.get("description", "")

    @classmethod
    def from_mat(cls, mat) -> "Model":
        """Build a Model from the struct loaded from a BiGG .mat file."""
        if isinstance(mat, Model):
            return mat

        if isinstance(mat, dict):
            get = mat.get
        elif mat_struct is not None and isinstance(mat, mat_struct):
            get = lambda name, default=None: getattr(mat, name, default)
        else:
            # A MATLAB struct may also arrive as a 1x1 object array.
            if isinstance(mat, np.ndarray) and mat.size == 1:
                return cls.from_mat(mat.ravel()[0])
            get = lambda name, default=None: getattr(mat, name, default)

        S = _as_matrix(get("S"))
        n_met, n_rxn = S.shape
        b = _as_vec(get("b"))
        c = _as_vec(get("c"))
        lb = _as_vec(get("lb"))
        ub = _as_vec(get("ub"))

        if b.size == 0:
            b = np.zeros(n_met)
        if c.size == 0:
            c = np.zeros(n_rxn)
        if lb.size == 0:
            lb = np.zeros(n_rxn)
        if ub.size == 0:
            ub = np.zeros(n_rxn)

        rxnGeneMat = get("rxnGeneMat")
        if rxnGeneMat is not None:
            rxnGeneMat = _as_matrix(rxnGeneMat)

        metCharge = get("metCharge")
        if metCharge is not None:
            metCharge = _as_vec(metCharge)

        rev = get("rev")
        if rev is not None:
            rev = _as_vec(rev)

        return cls(
            S=S,
            b=b,
            c=c,
            lb=lb,
            ub=ub,
            mets=_as_cellstr(get("mets")),
            metNames=_as_cellstr(get("metNames")),
            metFormulas=_as_cellstr(get("metFormulas")),
            metCharge=metCharge,
            rxns=_as_cellstr(get("rxns")),
            rxnNames=_as_cellstr(get("rxnNames")),
            subSystems=_as_cellstr(get("subSystems")),
            grRules=_as_cellstr(get("grRules")),
            genes=_as_cellstr(get("genes")),
            rev=rev,
            rxnGeneMat=rxnGeneMat,
            description=_as_str(get("description")),
        )

    @classmethod
    def load_mat(cls, filename: str, model_name: Optional[str] = None) -> "Model":
        """Load a model from a .mat file saved by MATLAB (BiGG format)."""
        data = sco.loadmat(filename, struct_as_record=False)
        keys = [key for key in data if not key.startswith("__")]
        if model_name is not None:
            if model_name not in data:
                raise KeyError(f"cannot find variable {model_name!r} in {filename}")
            key = model_name
        elif len(keys) == 1:
            key = keys[0]
        else:
            # Prefer a variable containing model fields.
            candidates = []
            for candidate in keys:
                value = data[candidate]
                if isinstance(value, np.ndarray) and value.size == 1:
                    value = value.ravel()[0]
                if mat_struct is not None and isinstance(value, mat_struct):
                    if "S" in getattr(value, "_fieldnames", []):
                        candidates.append(candidate)
            if len(candidates) != 1:
                raise ValueError(
                    "several variables found in the .mat file; "
                    "please specify the model name"
                )
            key = candidates[0]

        value = data[key]
        if isinstance(value, np.ndarray) and value.size == 1:
            value = value.ravel()[0]
        return cls.from_mat(value)

    def copy(self) -> "Model":
        return copy.deepcopy(self)

    def __repr__(self) -> str:  # pragma: no cover - convenience only
        return (f"Model({len(self.rxns)} reactions, {len(self.mets)} metabolites, "
                f"{len(self.genes)} genes)")
