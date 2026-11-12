"""Backward-compatible single-module entry point for :mod:`ratgenepy`.

The implementation has been split into the ``ratgenepy`` package.  This file
keeps older code working, for example::

    import RatGene as rg
    rg.RatGene(model, targetMet)

and also allows running the command line interface directly::

    python RatGene.py model.mat target_met
"""

from __future__ import annotations

import sys

from ratgenepy import *  # noqa: F401,F403
from ratgenepy import __all__ as _public_api
from ratgenepy import __version__  # noqa: F401
from ratgenepy._utils import (
    _as_cellstr,
    _as_matrix,
    _as_str,
    _as_vec,
    _blkdiag,
    _first_match,
    _gene_index,
    _normalize_vtype,
    _vstack,
)
from ratgenepy.candidates import _crossSet, _formBoolMat, _mergeSet  # noqa: F401
from ratgenepy.cli import _load_model_file, main
from ratgenepy.solvers import (  # noqa: F401
    available_solvers,
    cplex_solver,
    gurobi_solver,
    pool_solutions,
    scip_solver,
)
from ratgenepy.scip import (  # noqa: F401
    _SCIP_VARIABLES,
    _make_problem,
    _make_scip_variables,
    _maximize,
    _minimize,
    _pool_solutions,
    _solution_vector,
)

__all__ = list(_public_api)


if __name__ == "__main__":
    sys.exit(main())
