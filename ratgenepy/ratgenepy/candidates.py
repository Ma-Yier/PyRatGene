"""Candidate deletion/addition set construction and Boolean matrices."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ._utils import _gene_index

def _mergeSet(andSet, canDelCell, iterations, mode):
    if mode not in (0, 1):
        raise ValueError("mode only 1 or 0")

    if np.any(andSet == 1):
        iterations += 1
        while len(canDelCell) < iterations:
            canDelCell.append(None)
        canDelCell[iterations - 1] = np.asarray(andSet).reshape(1, -1)

    canDel = None
    for i in range(2, iterations + 1):
        if i == 2:
            set1 = canDelCell[0]
            set2 = canDelCell[1]
            canDel = _crossSet(set1, set2, mode)
        else:
            set1 = canDel
            set2 = canDelCell[i - 1]
            canDel = _crossSet(set1, set2, mode)

    if mode == 1 and iterations > 1:
        length = len(andSet)
        lengths = np.zeros((canDel.shape[0], 1))
        for i in range(canDel.shape[0]):
            lengths[i, 0] = np.sum(canDel[i, :])
            if length > lengths[i, 0]:
                length = lengths[i, 0]
        canDel = canDel[lengths[:, 0] == length, :]
    else:
        canDel = canDelCell[0]

    return canDel


def _crossSet(set1: np.ndarray, set2: np.ndarray, mode: int) -> np.ndarray:
    # The MATLAB implementation iterates over `size(set1,1)` and
    # `size(set2,1)` without a colon, i.e. it only combines the last rows.
    # This behaviour is kept to stay consistent with the reference code.
    length = set1.shape[1]
    finalSet = np.zeros((0, length))
    if mode == 0:
        for i in [set1.shape[0] - 1]:
            for j in [set2.shape[0] - 1]:
                a = set1[i, :] + set2[j, :]
                bb = np.zeros(length)
                bb[a > 0] = 1
                if np.sum(bb == 1) < length:
                    length = int(np.sum(bb == 1))
                    finalSet = bb.reshape(1, -1)
                elif np.sum(bb == 1) == length:
                    if finalSet.size == 0:
                        finalSet = bb.reshape(1, -1)
                    else:
                        c = finalSet @ bb
                        if not np.any(c == length):
                            finalSet = np.vstack((finalSet, bb))
                else:
                    continue
    elif mode == 1:
        lengthSet = np.zeros((0, 1))
        for i in [set1.shape[0] - 1]:
            for j in [set2.shape[0] - 1]:
                a = set1[i, :] + set2[j, :]
                bb = np.zeros(length)
                bb[a > 0] = 1
                lengthb = int(np.sum(bb == 1))
                if finalSet.size == 0:
                    cmpSet = np.zeros((0, length))
                else:
                    cmpSet = finalSet[lengthSet[:, 0] == lengthb, :]
                if i == 0 and j == 0:
                    finalSet = np.vstack((finalSet, bb))
                    lengthSet = np.vstack((lengthSet, lengthb))
                    continue
                if cmpSet.size == 0:
                    cmpSet = np.zeros((1, length))
                if np.any(cmpSet @ bb == lengthb):
                    continue
                finalSet = np.vstack((finalSet, bb))
                lengthSet = np.vstack((lengthSet, lengthb))
    return finalSet


def _formBoolMat(canDel: np.ndarray):
    numRela, numElem = canDel.shape
    matE = np.zeros((2 * numRela, numElem))
    matC = np.zeros((2 * numRela, numRela))
    matR = np.zeros((2 * numRela, 1))
    for i in range(numRela):
        matE[2 * i, canDel[i, :] == 1] = 1
        matE[2 * i + 1, canDel[i, :] == 1] = -1
        matC[2 * i, i] = -1
        matC[2 * i + 1, i] = np.sum(canDel[i, :])
        matR[2 * i, 0] = np.sum(canDel[i, :]) - 1
        matR[2 * i + 1, 0] = 0
    return matE, matC, matR


def getCandidateDel(parInfo, genes: Sequence[str], deletion: np.ndarray) -> np.ndarray:
    """Candidate gene deletion sets for one GPR clause."""
    parInfo = list(parInfo)
    nRel = len(parInfo)
    n_reactions = parInfo[-1][0]
    relation = parInfo[-1][1]
    parUnit = parInfo[-1][2]

    if relation == "and":
        canDelCell = []
        for gene in parUnit:
            if len(gene) >= 3 and gene[:3] == "aux":
                getID = [k for k, row in enumerate(parInfo) if row[0] == gene]
                if not getID:
                    raise ValueError(f"cannot locate {gene!r}")
                a = getCandidateDel(parInfo[:getID[0] + 1], genes, deletion)
                canDelCell.append(a)
            elif deletion[_gene_index(genes, gene)] == 0:
                a = np.zeros((1, len(deletion)))
                a[0, _gene_index(genes, gene)] = 1
                canDelCell.append(a)

        canDel = np.zeros((0, len(deletion)))
        length = len(deletion)
        for a in canDelCell:
            if np.sum(a[0, :] == 1) == 0:
                continue
            count = int(np.sum(a[0, :] == 1))
            if length > count:
                canDel = a
                length = count
            elif length == count:
                canDel = np.vstack((canDel, a))
        if canDel.shape[0] == 0:
            canDel = np.zeros((1, len(deletion)))
        return canDel

    if relation == "or":
        canDelCell = []
        iterations = 0
        andSet = np.zeros(len(deletion))
        sign = 0
        for gene in parUnit:
            if sign == 1:
                break
            if len(gene) >= 3 and gene[:3] == "aux":
                getID = [k for k, row in enumerate(parInfo) if row[0] == gene]
                if not getID:
                    raise ValueError(f"cannot locate {gene!r}")
                a = getCandidateDel(parInfo[:getID[0] + 1], genes, deletion)
                if np.sum(a[0, :] == 1) == 0:
                    andSet = np.zeros(len(deletion))
                    sign = 1
                else:
                    iterations += 1
                    if len(canDelCell) < iterations:
                        canDelCell.append(None)
                    canDelCell[iterations - 1] = a
            else:
                index = _gene_index(genes, gene)
                if deletion[index] == 0:
                    andSet[index] = 1
                elif deletion[index] == 1:
                    andSet = np.zeros(len(deletion))
                    sign = 1

        if sign == 1:
            return np.zeros((1, len(deletion)))
        if len(canDelCell) == 0 and not np.any(andSet == 1):
            return np.zeros((1, len(deletion)))
        if len(canDelCell) == 0:
            return andSet.reshape(1, -1)
        return _mergeSet(andSet, canDelCell, iterations, 1)

    return np.zeros((1, len(deletion)))


def getCandidateAdd(parInfo, genes: Sequence[str], addition: np.ndarray) -> np.ndarray:
    """Candidate gene addition sets for one GPR clause."""
    parInfo = list(parInfo)
    relation = parInfo[-1][1]
    parUnit = parInfo[-1][2]

    if relation == "or":
        canAddCell = []
        for gene in parUnit:
            if len(gene) >= 3 and gene[:3] == "aux":
                getID = [k for k, row in enumerate(parInfo) if row[0] == gene]
                if not getID:
                    raise ValueError(f"cannot locate {gene!r}")
                a = getCandidateAdd(parInfo[:getID[0] + 1], genes, addition)
                canAddCell.append(a)
            elif addition[_gene_index(genes, gene)] == 1:
                a = np.zeros((1, len(addition)))
                a[0, _gene_index(genes, gene)] = 1
                canAddCell.append(a)

        canAdd = np.zeros((0, len(addition)))
        length = len(addition)
        for a in canAddCell:
            if np.sum(a[0, :] == 1) == 0:
                continue
            count = int(np.sum(a[0, :] == 1))
            if length > count:
                canAdd = a
                length = count
            elif length == count:
                canAdd = np.vstack((canAdd, a))
        if canAdd.shape[0] == 0:
            canAdd = np.zeros((1, len(addition)))
        return canAdd

    if relation == "and":
        canAddCell = []
        iterations = 0
        andSet = np.zeros(len(addition))
        sign = 0
        for gene in parUnit:
            if sign == 1:
                break
            if len(gene) >= 3 and gene[:3] == "aux":
                getID = [k for k, row in enumerate(parInfo) if row[0] == gene]
                if not getID:
                    raise ValueError(f"cannot locate {gene!r}")
                a = getCandidateAdd(parInfo[:getID[0] + 1], genes, addition)
                if np.sum(a[0, :] == 1) == 0:
                    andSet = np.zeros(len(addition))
                    sign = 1
                else:
                    iterations += 1
                    if len(canAddCell) < iterations:
                        canAddCell.append(None)
                    canAddCell[iterations - 1] = a
            else:
                index = _gene_index(genes, gene)
                if addition[index] == 1:
                    andSet[index] = 1
                elif addition[index] == 0:
                    andSet = np.zeros(len(addition))
                    sign = 1

        if sign == 1:
            return np.zeros((1, len(addition)))
        if len(canAddCell) == 0 and not np.any(andSet == 1):
            return np.zeros((1, len(addition)))
        if len(canAddCell) == 0:
            return andSet.reshape(1, -1)
        return _mergeSet(andSet, canAddCell, iterations, 1)

    return np.zeros((1, len(addition)))
