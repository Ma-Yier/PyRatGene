"""Ratio-based gene and reaction modification algorithms."""

from __future__ import annotations

import time

import numpy as np

from ._resources import parallel_limits
from .gpr import constructMatrix
from .model import Model
from .solvers import pool_solutions, solver
from .verification import verifyGeneKnock, verifyKnock

def ratioGene(model: Model, id_biomass: int, id_target: int, TMGR: float,
              max_loop: int, gap: float, numMultiStrat: int,
              timeLimit: float, backend: str = "scip",
              integerTolerance=None, threads: int = 1,
              memory: float = float("inf")):
    """Gene deletion or deletion/addition strategy search (ratioGene.m).

    ``threads`` limits the number of solver threads (default ``1``) and
    ``memory`` limits the total solver memory in MB (default ``inf``, meaning
    no limit).
    """
    (lessMatrixLeft, lessMatrixRight, equalMatrixLeft, equalMatrixRight,
     nGpr, nAux, nGen, indGPR) = constructMatrix(model)

    alpha = 0.0
    n = model.S.shape[1]
    row_ratio = equalMatrixLeft.shape[0]

    A = lessMatrixLeft
    b = lessMatrixRight

    ratioCons = np.zeros((1, n + nGpr + nAux + nGen))
    ratioCons[0, id_target] = 1
    ratioCons[0, id_biomass] = 0

    Aeq = np.vstack((equalMatrixLeft, ratioCons))
    beq = np.vstack((equalMatrixRight, np.zeros((1, 1))))

    lb = np.concatenate((model.lb, np.zeros(nGpr + nAux + nGen)))
    ub = np.concatenate((model.ub, np.ones(nGpr + nAux + nGen)))

    f = np.zeros(n + nGpr + nAux + nGen)
    f[n:n + nGpr] = TMGR
    f[id_biomass] = -1

    vtype = "C" * n + "B" * (nGpr + nAux + nGen)

    x_target = 0.0
    knockout = np.array([])
    start = time.time()

    for _ in range(max_loop):
        alpha += gap
        Aeq[row_ratio, id_biomass] = -alpha
        Aeq[row_ratio, id_target] = 1

        remaining = timeLimit - (time.time() - start)
        if remaining <= 0:
            break
        milpTime = 100.0 if not np.isfinite(timeLimit) else max(1.0, min(100.0, remaining))

        handle, best = solver(
            f, A=A, b=b, Aeq=Aeq, beq=beq, lb=lb, ub=ub,
            opt="milp", vtype=vtype, timeLimit=milpTime,
            pool=numMultiStrat, backend=backend,
            integerTolerance=integerTolerance,
            threads=threads, memory=memory,
        )

        if best is not None:
            # The variables are ordered as [flux; reaction indicators;
            # auxiliary indicators; genes].
            gene_offset = n + nGpr + nAux
            pool = pool_solutions(handle, numMultiStrat)
            if len(pool) == 0:
                pool = [best]
            for solution in pool:
                if solution.size != n + nGpr + nAux + nGen:
                    x_target = 0.0
                    continue
                gene_values = solution[gene_offset:gene_offset + nGen]
                x_target, knockout_gene = verifyGeneKnock(
                    model, gene_values, indGPR, id_biomass, id_target,
                    model.lb[id_biomass], backend=backend,
                    threads=threads, memory=memory,
                )
                if x_target > 0:
                    return x_target, alpha, knockout_gene
                x_target = 0.0

        if (time.time() - start) > timeLimit:
            break

    return x_target, alpha, knockout


def ratioGenePar(model: Model, id_biomass: int, id_target: int, TMGR: float,
                 max_loop: int, gap: float, numMultiStrat: int,
                 timeLimit: float, backend: str = "scip",
                 integerTolerance=None, threads: int = 1,
                 memory: float = float("inf")):
    """Parallel-resource-safe version of :func:`ratioGene`.

    The signature is identical to :func:`ratioGene`.  Before running the
    algorithm, common BLAS/OpenMP/NumExpr environment variables are pinned to
    ``threads`` so that NumPy/SciPy matrix operations do not use all available
    cores.  When ``memory`` is finite, an operating-system-level limit is
    applied for the duration of the call: ``resource.RLIMIT_AS`` on POSIX and
    a Job Object process-memory limit on Windows.  Solver-level limits are
    also passed to the selected backend.

    ``memory`` is expressed in megabytes; the default ``inf`` disables the
    hard limit.
    """
    with parallel_limits(threads=threads, memory=memory):
        return ratioGene(
            model, id_biomass, id_target, TMGR, max_loop, gap,
            numMultiStrat, timeLimit, backend=backend,
            integerTolerance=integerTolerance, threads=threads,
            memory=memory,
        )


def ratioMethod(model: Model, id_biomass: int, id_target: int,
                max_loop: int, gap: float, timeLimit: float,
                backend: str = "scip", threads: int = 1,
                memory: float = float("inf")):
    """Reaction deletion or deletion/addition strategy search (ratioMethod.m)."""
    alpha = 0.0
    m, n = model.S.shape
    knockout = np.array([])
    x_target = 0.0
    start = time.time()

    for _ in range(max_loop):
        alpha += gap
        ratio = np.zeros(n)
        ratio[id_biomass] = -alpha
        ratio[id_target] = 1

        Aeq = np.hstack((np.vstack((model.S, ratio)),
                         np.zeros((m + 1, n))))
        beq = np.append(model.b, 0.0)

        A = np.vstack((np.hstack((np.eye(n), -np.eye(n))),
                       np.hstack((-np.eye(n), -np.eye(n)))))
        b = np.zeros(2 * n)
        c = np.concatenate((np.zeros(n), np.ones(n)))
        lb = np.concatenate((model.lb, np.zeros(n)))
        ub = np.concatenate((model.ub, 999999.0 * np.ones(n)))

        remaining = timeLimit - (time.time() - start)
        if remaining <= 0:
            break

        _, x = solver(c, A=A, b=b, Aeq=Aeq, beq=beq, lb=lb, ub=ub,
                      opt="lp", timeLimit=remaining, backend=backend,
                      threads=threads, memory=memory)

        if x is not None and x.size == 2 * n:
            x_target, tknockout = verifyKnock(
                model, x[:n], id_biomass, id_target, model.lb[id_biomass],
                backend=backend, threads=threads, memory=memory
            )
            if x_target > 0:
                knockout = tknockout
                break
        else:
            x_target = 0.0

        if (time.time() - start) > timeLimit:
            break

    return x_target, alpha, knockout


def RatMethodRxn(model: Model, biomassID: int, targetID: int, maxLoop: int,
                 gap: float, timeLimit: float, backend: str = "scip",
                 threads: int = 1, memory: float = float("inf")):
    """Compatibility wrapper for the reaction-level ratio method.

    Returns ``(x_target, knockout, alpha)`` as in the original Python draft.
    """
    x_target, alpha, knockout = ratioMethod(
        model, biomassID, targetID, maxLoop, gap, timeLimit,
        backend=backend, threads=threads, memory=memory
    )
    return x_target, knockout, alpha


def verifyRxn(model: Model, knot, biomassID: int, targeID: int,
              LBbiomass: float, backend: str = "scip", threads: int = 1,
              memory: float = float("inf")):
    """Compatibility wrapper around :func:`verifyKnock`.

    Returns ``(x_target, deleted_reactions)``.
    """
    x_target, knockout = verifyKnock(model, knot, biomassID, targeID,
                                     LBbiomass, backend=backend,
                                     threads=threads, memory=memory)
    return x_target, np.where(np.asarray(knockout) == 0)[0]
