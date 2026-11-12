"""GPR parsing, Boolean-matrix construction and GPR evaluation."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ._utils import _as_str, _as_vec, _gene_index
from .model import Model

def maxParseGPR(model: Model):
    """Parse GPR rules into linear Boolean relationships.

    This follows ``maxParseGPR.m``: every Boolean clause (one parenthesis) is
    represented by two inequalities, therefore the number of auxiliary
    variables is not minimised.
    """
    nRxn = len(model.rxns)
    nGen = len(model.genes)
    nAux = 1
    parInfo = [None] * nRxn
    nRelation = 0
    nEqual = 0
    gprLabel = np.ones(nRxn, dtype=float)

    for i in range(nRxn):
        gpr = _as_str(model.grRules[i]) if i < len(model.grRules) else ""
        if gpr:
            if " " in gpr:
                gprStr = "(" + gpr + ")"
                maxUnit = gprStr.count("(")
                nRelation += 2 * maxUnit
                parInfoUnit = [[None, None, None] for _ in range(maxUnit)]
                unitLabel = 0

                while unitLabel < maxUnit:
                    rpoint = gprStr.find(")")
                    if rpoint < 0:
                        raise ValueError(f"cannot parse GPR rule: {gpr!r}")
                    lpoint = gprStr.rfind("(", 0, rpoint)
                    if lpoint < 0:
                        raise ValueError(f"cannot parse GPR rule: {gpr!r}")
                    gprUnit = gprStr[lpoint + 1:rpoint]

                    if unitLabel == maxUnit - 1:
                        name = "real" + str(i + 1)
                    else:
                        name = "aux" + str(nAux)
                        nAux += 1

                    gprStr = gprStr[:lpoint] + name + gprStr[rpoint + 1:]
                    elements = [element for element in gprUnit.split(" ") if element]
                    if len(elements) < 3 or elements[1] not in ("and", "or"):
                        raise ValueError(f"cannot parse GPR rule: {gpr!r}")
                    parInfoUnit[unitLabel] = [name, elements[1],
                                              elements[0::2]]
                    unitLabel += 1
            else:
                parInfoUnit = [["real" + str(i + 1), "", [gpr]]]
                nEqual += 1
        else:
            parInfoUnit = [[None, None, None]]
            gprLabel[i] = 0
        parInfo[i] = parInfoUnit

    nAux = nAux - 1
    return parInfo, nRxn, nGen, nAux, nRelation, nEqual, gprLabel


def locateGene(nRxn: int, nGpr: int, nAux: int, indRea: np.ndarray,
               objGene: str, model: Model) -> int:
    """Find the 0-based column of a gene/clause in the MILP variable vector."""
    geneNameLen = len(objGene)

    if geneNameLen > 2:
        head = objGene[:3]
        if head == "rea":
            ind = int(objGene[4:geneNameLen]) - 1
            return int(indRea[ind]) + nRxn
        if head == "aux":
            return int(objGene[3:geneNameLen]) - 1 + nRxn + nGpr

    return _gene_index(model.genes, objGene) + nRxn + nGpr + nAux


def constructMatrix(model: Model):
    """Construct the equality/inequality matrices describing GPR logic."""
    (parInfo, nRxn, nGen, nAux, nRelation, nEqual,
     gprLabel) = maxParseGPR(model)

    indGPR = np.where(gprLabel != 0)[0]
    nGpr = len(indGPR)
    indRea = np.zeros(nRxn, dtype=int)
    for i in range(nGpr):
        indRea[indGPR[i]] = i

    nMet = model.S.shape[0]
    nVar = nRxn + nGpr + nGen + nAux

    gprMatrixLeft = np.zeros((nRelation, nVar))
    gprMatrixRight = np.zeros((nRelation, 1))

    gprLabelDiag = np.diag(gprLabel)
    grCorrelationRxn = gprLabelDiag[indGPR, :]
    boundMatrixLeft = np.vstack((
        np.hstack((grCorrelationRxn,
                   -np.diag(model.ub[indGPR]),
                   np.zeros((nGpr, nGen + nAux)))),
        np.hstack((-grCorrelationRxn,
                   np.diag(model.lb[indGPR]),
                   np.zeros((nGpr, nGen + nAux)))),
    ))
    boundMatrixRight = np.zeros((2 * nGpr, 1))

    gprEqualLeft = np.zeros((nEqual, nVar))
    gprEqualRight = np.zeros((nEqual, 1))

    fluxEqualLeft = np.hstack((
        model.S, np.zeros((nMet, nGpr + nGen + nAux))
    ))
    fluxEqualRight = np.zeros((nMet, 1))

    rowPoint1 = 0
    rowPoint2 = 0

    for i in range(len(parInfo)):
        parInfoUnit = parInfo[i]

        if parInfoUnit[0][0] is None:
            continue

        maxUnit = len(parInfoUnit)

        if maxUnit == 1 and not parInfoUnit[0][1]:
            geneID = _gene_index(model.genes, parInfoUnit[0][2][0])
            gprEqualLeft[rowPoint1, nRxn + indRea[i]] = 1
            gprEqualLeft[rowPoint1, nRxn + nGpr + nAux + geneID] = -1
            rowPoint1 += 1
        else:
            for j in range(maxUnit):
                name, logicRelation, ctrlGenes = parInfoUnit[j]
                n = len(ctrlGenes)
                objLoc = locateGene(nRxn, nGpr, nAux, indRea, name, model)

                if logicRelation == "and":
                    gprMatrixLeft[rowPoint2, objLoc] = -1
                    gprMatrixLeft[rowPoint2 + 1, objLoc] = n
                    gprMatrixRight[rowPoint2, 0] = n - 1
                elif logicRelation == "or":
                    gprMatrixLeft[rowPoint2, objLoc] = -n
                    gprMatrixLeft[rowPoint2 + 1, objLoc] = 1
                else:
                    raise ValueError("parsing GPR error")

                for gene in ctrlGenes:
                    geneLoc = locateGene(nRxn, nGpr, nAux, indRea, gene, model)
                    gprMatrixLeft[rowPoint2, geneLoc] = 1
                    gprMatrixLeft[rowPoint2 + 1, geneLoc] = -1
                rowPoint2 += 2

    lessMatrixLeft = np.vstack((boundMatrixLeft, gprMatrixLeft))
    lessMatrixRight = np.vstack((boundMatrixRight, gprMatrixRight))
    equalMatrixLeft = np.vstack((fluxEqualLeft, gprEqualLeft))
    equalMatrixRight = np.vstack((fluxEqualRight, gprEqualRight))

    return (lessMatrixLeft, lessMatrixRight,
            equalMatrixLeft, equalMatrixRight,
            nGpr, nAux, nGen, indGPR)


def verifyRatioGene(model: Model, geneValue) -> np.ndarray:
    """Evaluate GPR rules for a 0/1 gene modification strategy."""
    geneValue = _as_vec(geneValue)
    grRules = [str(rule).replace("or", "+").replace("and", "*")
               for rule in model.grRules]
    genes = list(model.genes)
    num_rxn = len(grRules)
    geneKnock = np.zeros(num_rxn)

    order = sorted(range(len(genes)), key=lambda k: genes[k], reverse=True)
    for k in order:
        value = geneValue[k] if k < len(geneValue) else 1
        value_str = str(int(round(value)))
        grRules = [rule.replace(genes[k], value_str) for rule in grRules]

    for j, rule in enumerate(grRules):
        if rule.strip() == "":
            continue
        try:
            value = eval(rule, {"__builtins__": {}}, {})
        except Exception as exc:
            raise ValueError(f"cannot evaluate GPR rule {model.grRules[j]!r}") from exc
        if value < 0.9:
            geneKnock[j] = 1

    return geneKnock
