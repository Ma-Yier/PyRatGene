"""Two-step reduction of the size of a modification strategy."""

from __future__ import annotations

import numpy as np

from ._utils import _as_str, _as_vec, _blkdiag, _vstack
from .candidates import _formBoolMat, getCandidateAdd, getCandidateDel
from .gpr import maxParseGPR, verifyRatioGene
from .model import Model
from .solvers import solver

def rmRedundancy(model: Model, knockouts, id_target: int, id_biomass: int,
                 split, backend: str = "scip"):
    """Two-step approach reducing the size of a modification strategy.

    ``split`` is ``[0, 0]`` for a pure deletion problem or
    ``[number_of_core_reactions, number_of_core_genes]`` for an
    addition problem.
    """
    afterNumKnockouts = -1
    preNumKnockouts = -1
    newKnockouts = np.array([np.inf])

    knockouts = _as_vec(knockouts).copy()
    split = _as_vec(split)

    if knockouts.size == len(model.genes):
        pass
    elif knockouts.size == len(model.rxns):
        print("rxn knockouts")
        return preNumKnockouts, afterNumKnockouts, newKnockouts
    else:
        print("neither rxn knockout nor gene knockout")
        return preNumKnockouts, afterNumKnockouts, newKnockouts

    if split[0] > 0:
        model1 = model.copy()
        model2 = model.copy()
        model1.genes = list(model.genes[:int(split[1])])
        model1.grRules = list(model.grRules[:int(split[0])])
        model1.rxns = list(model.rxns[:int(split[0])])
        model2.genes = list(model.genes[int(split[1]):])
        model2.grRules = list(model.grRules[int(split[0]):])
        model2.rxns = list(model.rxns[int(split[0]):])
        endRxn = int(split[0])
        endGene = int(split[1])
        additions = knockouts[int(split[1]):].copy()
        knockouts = knockouts[:int(split[1])].copy()
    else:
        model1 = model.copy()
        endRxn = len(model1.rxns)
        endGene = len(model1.genes)
        additions = np.array([])

    parInfo, _, _, _, _, _, _ = maxParseGPR(model1)
    knockouts[knockouts > 0.1] = 1
    var_x = verifyRatioGene(model1, knockouts)

    newKnockouts = _reduce_one(
        model1, knockouts, var_x, endRxn, endGene, parInfo,
        deletion=True, backend=backend
    )
    preDel = int(np.sum(knockouts == 0))
    afterDel = int(np.sum(newKnockouts == 0))

    if split[0] == 0:
        return preDel, afterDel, newKnockouts

    parInfo2, _, _, _, _, _, _ = maxParseGPR(model2)
    additions[additions > 0.1] = 1
    var_a = verifyRatioGene(model2, additions)

    newAdditions = _reduce_one(
        model2, additions, var_a, len(model2.rxns), len(model2.genes),
        parInfo2, deletion=False, backend=backend
    )
    preAdd = int(np.sum(additions == 1))
    afterAdd = int(np.sum(newAdditions == 1))

    return (np.array([preDel, preAdd]), np.array([afterDel, afterAdd]),
            np.concatenate((newKnockouts, newAdditions)))


def _reduce_one(model: Model, knockouts: np.ndarray, var_x: np.ndarray,
                endRxn: int, endGene: int, parInfo, deletion: bool,
                backend: str = "scip"):
    genes = list(model.genes[:endGene])
    grRules = model.grRules

    rxnNum = 0
    b = np.zeros((0, 1))
    splitRxn = np.zeros((0, 0))
    splitGene = np.zeros((0, 0))
    splitAux = np.zeros((0, 0))
    equalRxn = np.zeros((0, 0))
    equalGene = np.zeros((0, 0))

    for i in range(endRxn):
        rule = _as_str(grRules[i])
        if not rule:
            continue
        # Deletion strategies consider reactions that are knocked out
        # (var_x == 1); addition strategies consider reactions that become
        # active (var_x == 0).  See the two loops in rmRedundancy.m.
        if deletion:
            if var_x[i] == 0:
                continue
        elif var_x[i] == 1:
            continue
        rxnNum += 1

        if rule in model.genes:
            equalRxn = _blkdiag(equalRxn, np.ones((1, 1)))
            eg = np.zeros((1, len(knockouts)))
            eg[0, model.genes.index(rule)] = -1
            equalGene = _vstack(equalGene, eg)
            continue

        if deletion:
            canDel = getCandidateDel(parInfo[i], genes, knockouts)
        else:
            canDel = getCandidateAdd(parInfo[i], genes, knockouts)
        if np.sum(canDel) == 0:
            raise ValueError("error del" if deletion else "error add")

        if canDel.shape[0] == 1:
            interMatGene = np.vstack((canDel, -canDel))
            interMatRxn = np.array([[-1.0], [np.sum(canDel)]])
            interMatb = np.array([[np.sum(canDel) - 1.0], [0.0]])

            b = _vstack(b, interMatb)
            splitGene = _vstack(splitGene, interMatGene)
            splitAux = _vstack(splitAux, np.zeros((2, splitAux.shape[1])))
            splitRxn = _blkdiag(splitRxn, interMatRxn)
        else:
            matGene, matAux, matb = _formBoolMat(canDel)
            nAux = matAux.shape[1]

            interMatGene = np.vstack((np.zeros((2, matGene.shape[1])), matGene))
            interMatb = np.vstack((np.zeros((2, 1)), matb))
            interMatAux = np.vstack((-np.ones((1, nAux)),
                                     np.ones((1, nAux)), matAux))
            interMatRxn = np.vstack((np.ones((1, 1)),
                                     -float(nAux) * np.ones((1, 1)),
                                     np.zeros((matGene.shape[0], 1))))

            b = _vstack(b, interMatb)
            splitGene = _vstack(splitGene, interMatGene)
            splitAux = _blkdiag(splitAux, interMatAux)
            splitRxn = _blkdiag(splitRxn, interMatRxn)

    if deletion:
        mask = knockouts == 0
    else:
        mask = knockouts == 1

    splitGene2 = splitGene[:, mask] if splitGene.size else splitGene
    equalGene2 = equalGene[:, mask] if equalGene.size else equalGene
    a = splitRxn.shape[1] + equalRxn.shape[1]
    bb = splitAux.shape[1]
    cnum = splitGene2.shape[1]

    if cnum == 0:
        return knockouts.copy()

    Aeq = np.hstack((equalRxn,
                     np.zeros((equalRxn.shape[0], splitRxn.shape[1] + bb)),
                     equalGene2))
    beq = np.zeros(equalRxn.shape[0])
    A = np.hstack((np.zeros((splitRxn.shape[0], equalRxn.shape[1])),
                   splitRxn, splitAux, splitGene2))
    f = np.concatenate((np.zeros(a + bb), np.ones(cnum)))
    lb = np.concatenate((np.ones(a), np.zeros(bb + cnum)))
    ub = np.ones(a + bb + cnum)
    vtype = "B" * (a + bb + cnum)

    _, x = solver(f, A=A, b=b.ravel(), Aeq=Aeq, beq=beq, lb=lb, ub=ub,
                  opt="milp", vtype=vtype, backend=backend)

    if x is None or cnum == 0:
        return knockouts.copy()

    # Round the binary selection variables; SCIP may return values such as
    # 0.9999999 depending on the integrality tolerance.
    newSel = np.rint(x[a + bb:a + bb + cnum]).astype(float)
    newSel = newSel - 1
    newSel[newSel == -1] = 1

    if deletion:
        # rmRedundancy.m: del = knockouts - 1 (with -1 mapped to 1), then
        # newKnockouts = knockouts + del * newDel.
        base = knockouts - 1
        base[base == -1] = 1
        base = np.diag(base)[:, mask]
        newStrategy = knockouts + base @ newSel
    else:
        # rmRedundancy.m: add = diag(additions), then
        # newAdditions = additions - add * newAdd.
        base = np.diag(knockouts)[:, mask]
        newStrategy = knockouts - base @ newSel

    var_x2 = verifyRatioGene(model, newStrategy)
    index1 = np.where(var_x == 1)[0]
    index2 = np.where(var_x2 == 1)[0]
    if not np.array_equal(index1, index2):
        newStrategy = knockouts.copy()

    return newStrategy
