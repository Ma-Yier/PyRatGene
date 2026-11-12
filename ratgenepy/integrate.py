"""Integration of a core model and an addition model."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ._utils import _as_matrix, _as_str, _gene_index
from .model import Model

def findDiffGene(model: Model, diff_rxns: Sequence[int]) -> np.ndarray:
    """Find all genes that belong to a set of reactions."""
    genes = []
    for rxn in diff_rxns:
        grRule = _as_str(model.grRules[int(rxn)])
        if not grRule:
            continue
        if " " in grRule:
            grRule = grRule.replace("(", "").replace(")", "")
            for item in grRule.split():
                if item in ("and", "or"):
                    continue
                genes.append(item)
        else:
            genes.append(grRule)

    genesid = []
    point = 0
    while point < len(genes):
        pattern = genes[point]
        genesid.append(_gene_index(model.genes, pattern))
        genes = [gene for k, gene in enumerate(genes) if not (gene == pattern and k != point)]
        point += 1

    genesid = np.sort(np.asarray(genesid, dtype=int))
    return genesid


def uniteGeneModel(model1: Model, model2: Model) -> Model:
    """Unite two models into one, used for gene deletion/addition problems."""
    model1 = model1.copy() if isinstance(model1, Model) else Model.from_mat(model1)
    model2 = model2.copy() if isinstance(model2, Model) else Model.from_mat(model2)

    m, n = model2.S.shape
    a = len(model2.genes)
    x, y = model1.S.shape
    b = len(model1.genes)

    metPair = -np.ones((m, 2), dtype=int)
    rxnPair = -np.ones((n, 2), dtype=int)
    genePair = -np.ones((a, 2), dtype=int)

    for i in range(m):
        for j in range(x):
            if str(model2.metNames[i]) == str(model1.metNames[j]):
                metPair[i, :] = [i, j]
                break

    for i in range(n):
        for j in range(y):
            if str(model2.rxnNames[i]) == str(model1.rxnNames[j]):
                rxnPair[i, :] = [i, j]
                break

    for i in range(a):
        for j in range(b):
            if str(model2.genes[i]) == str(model1.genes[j]):
                genePair[i, :] = [i, j]
                break

    diff_mets = np.where(metPair[:, 0] == -1)[0]
    diff_rxns = np.where(rxnPair[:, 0] == -1)[0]
    num_diff_mets = len(diff_mets)
    num_diff_rxns = len(diff_rxns)
    same_mets = np.where(metPair[:, 0] != -1)[0]
    same_rxns = np.where(rxnPair[:, 0] != -1)[0]
    num_same_mets = len(same_mets)
    num_same_rxns = len(same_rxns)

    diff_genes = findDiffGene(model2, diff_rxns)
    num_diff_genes = len(diff_genes)
    same_genes = np.where(genePair[:, 0] != -1)[0]
    num_same_genes = len(same_genes)

    universal = model1.copy()
    diff_matrix = model2.S[:, diff_rxns] if num_diff_rxns else np.zeros((m, 0))

    rightMatrix = np.zeros((x, num_diff_rxns))
    for idx in same_mets:
        rightMatrix[metPair[idx, 1], :] = diff_matrix[metPair[idx, 0], :]

    underMatrix = np.zeros((num_diff_mets, y))
    rightUnderMatrix = model2.S[np.ix_(diff_mets, diff_rxns)]

    universal.S = np.block([[model1.S, rightMatrix],
                            [underMatrix, rightUnderMatrix]])

    universal.lb = np.concatenate((model1.lb, model2.lb[diff_rxns]))
    universal.ub = np.concatenate((model1.ub, model2.ub[diff_rxns]))
    universal.b = np.concatenate((model1.b, np.zeros(num_diff_mets)))
    universal.c = np.concatenate((model1.c, np.zeros(num_diff_rxns)))
    if model1.rev is not None and model2.rev is not None:
        universal.rev = np.concatenate((model1.rev, model2.rev[diff_rxns]))
    else:
        universal.rev = None

    universal.mets = list(model1.mets) + [model2.mets[i] for i in diff_mets]
    universal.metNames = list(model1.metNames) + [model2.metNames[i] for i in diff_mets]
    universal.metFormulas = (list(model1.metFormulas)
                             + [model2.metFormulas[i] for i in diff_mets])
    if model1.metCharge is not None and model2.metCharge is not None:
        universal.metCharge = np.concatenate(
            (model1.metCharge, model2.metCharge[diff_mets])
        )
    else:
        universal.metCharge = None

    universal.grRules = list(model1.grRules) + [model2.grRules[i] for i in diff_rxns]
    universal.rxns = list(model1.rxns) + [model2.rxns[i] for i in diff_rxns]
    universal.rxnNames = list(model1.rxnNames) + [model2.rxnNames[i] for i in diff_rxns]
    universal.subSystems = (list(model1.subSystems)
                            + [model2.subSystems[i] for i in diff_rxns])

    universal.genes = list(model1.genes) + [model2.genes[i] for i in diff_genes]

    # rxnGeneMat
    if model1.rxnGeneMat is None:
        mat1 = np.zeros((y, b))
    else:
        mat1 = _as_matrix(model1.rxnGeneMat)
    if model2.rxnGeneMat is None:
        mat2 = np.zeros((n, a))
    else:
        mat2 = _as_matrix(model2.rxnGeneMat)

    southMatrix = np.zeros((num_diff_rxns, b))
    diffRxnMat = mat2[diff_rxns, :]
    for idx in same_genes:
        southMatrix[:, genePair[idx, 1]] = diffRxnMat[:, genePair[idx, 0]]

    eastMatrix = np.zeros((y, num_diff_genes))
    diffGeneMat = mat2[:, diff_genes]
    for idx in same_rxns:
        eastMatrix[rxnPair[idx, 1], :] = diffGeneMat[rxnPair[idx, 0], :]

    southEastMatrix = mat2[np.ix_(diff_rxns, diff_genes)]
    universal.rxnGeneMat = np.block([[mat1, eastMatrix],
                                     [southMatrix, southEastMatrix]])

    universal.description = (f"United models. Core: {model1.description}, "
                             f"Addition: {model2.description}")
    return universal
