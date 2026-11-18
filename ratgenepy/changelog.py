"""Human- and machine-readable change logs for RatGene modifications.

The log records every reaction whose bounds are changed while a RatGene
strategy is constructed.  Each entry contains the reaction identifier, the
reaction formula, its compartment, the old and new bounds, and the reason for
the change.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np

# Common BiGG compartment abbreviations.  Unknown abbreviations are kept as-is.
COMPARTMENT_NAMES = {
    "c": "cytosol",
    "e": "extracellular",
    "m": "mitochondrion",
    "n": "nucleus",
    "r": "endoplasmic reticulum",
    "g": "Golgi apparatus",
    "l": "lysosome",
    "p": "peroxisome",
    "x": "peroxisome",
    "v": "vacuole",
    "s": "secreted",
    "u": "unknown",
}


class ChangeLogEntry(dict):
    """A single bound change.

    Entries behave like ordinary dictionaries (so they can be serialized
    directly) and also expose their keys as attributes, e.g.
    ``entry.reaction_id``.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:  # pragma: no cover - attribute protocol
            raise AttributeError(name) from exc


class ChangeLog(list):
    """An ordered collection of :class:`ChangeLogEntry` dictionaries.

    ``ChangeLog`` is a ``list`` subclass, therefore it can be iterated,
    indexed and passed to code that expects a list of dictionaries.  Extra
    information about the run (for example which target outlet was selected)
    is available in :attr:`metadata`.
    """

    def __init__(self, entries: Optional[Iterable[Mapping]] = None):
        super().__init__()
        self.metadata = {}
        if entries is not None:
            for entry in entries:
                self.append(entry if isinstance(entry, ChangeLogEntry) else ChangeLogEntry(entry))

    # ------------------------------------------------------------------
    # recording helpers
    # ------------------------------------------------------------------
    def record(self, reaction_id: str, formula: str, compartment: str,
               old_bounds: Sequence[float], new_bounds: Sequence[float],
               reason: str, *, reaction_name: str = "",
               change_type: str = "bounds", role: Optional[str] = None,
               compartment_name: Optional[str] = None,
               introduced: bool = False,
               genes: Optional[Sequence[str]] = None) -> ChangeLogEntry:
        """Append a fully specified entry."""
        old_lb, old_ub = _bound_pair(old_bounds)
        new_lb, new_ub = _bound_pair(new_bounds)
        entry = ChangeLogEntry(
            reaction_id=str(reaction_id),
            reaction_name=str(reaction_name or ""),
            formula=str(formula or ""),
            compartment=str(compartment or ""),
            compartment_name=str(
                compartment_name
                if compartment_name is not None
                else compartment_names(compartment)
            ),
            old_bounds=(old_lb, old_ub),
            new_bounds=(new_lb, new_ub),
            old_lb=old_lb,
            old_ub=old_ub,
            new_lb=new_lb,
            new_ub=new_ub,
            reason=str(reason),
            change_type=str(change_type),
            role=role,
            introduced=bool(introduced),
            genes=list(genes or []),
        )
        self.append(entry)
        return entry

    def add_bounds_change(self, model, reaction_index: int,
                          old_bounds: Sequence[float],
                          new_bounds: Sequence[float], reason: str, *,
                          change_type: str = "bounds",
                          role: Optional[str] = None,
                          genes: Optional[Sequence[str]] = None,
                          introduced: bool = False,
                          force: bool = False) -> Optional[ChangeLogEntry]:
        """Record a change for a reaction of ``model``.

        Unchanged bounds are ignored unless ``force`` is true.
        """
        reaction_index = int(reaction_index)
        old_lb, old_ub = _bound_pair(old_bounds)
        new_lb, new_ub = _bound_pair(new_bounds)
        if not force and old_lb == new_lb and old_ub == new_ub:
            return None
        return self.record(
            reaction_id=_reaction_id(model, reaction_index),
            reaction_name=_reaction_name(model, reaction_index),
            formula=reaction_formula(model, reaction_index),
            compartment=reaction_compartment(model, reaction_index),
            old_bounds=(old_lb, old_ub),
            new_bounds=(new_lb, new_ub),
            reason=reason,
            change_type=change_type,
            role=role,
            genes=genes,
            introduced=introduced,
        )

    def add_introduced_reaction(self, model, reaction_index: int, reason: str, *,
                                role: Optional[str] = None,
                                change_type: str = "reaction_introduced",
                                old_bounds: Sequence[float] = (0.0, 0.0)):
        """Record a reaction that did not exist before this run."""
        reaction_index = int(reaction_index)
        return self.record(
            reaction_id=_reaction_id(model, reaction_index),
            reaction_name=_reaction_name(model, reaction_index),
            formula=reaction_formula(model, reaction_index),
            compartment=reaction_compartment(model, reaction_index),
            old_bounds=old_bounds,
            new_bounds=(model.lb[reaction_index], model.ub[reaction_index]),
            reason=reason,
            change_type=change_type,
            role=role,
            introduced=True,
        )

    # ------------------------------------------------------------------
    # serialization
    # ------------------------------------------------------------------
    def to_dicts(self):
        """Return the entries as plain dictionaries."""
        return [dict(entry) for entry in self]

    def to_json(self, path: Optional[str] = None, indent: int = 2):
        """Return JSON text and optionally write it to ``path``."""
        text = json.dumps(self.to_dicts(), indent=indent, ensure_ascii=False)
        if path is not None:
            Path(path).write_text(text, encoding="utf-8")
        return text

    def to_csv(self, path: Optional[str] = None):
        """Return CSV text and optionally write it to ``path``."""
        import io

        buffer = io.StringIO()
        fieldnames = [
            "reaction_id", "reaction_name", "formula", "compartment",
            "compartment_name", "old_lb", "old_ub", "new_lb", "new_ub",
            "reason", "change_type", "role", "introduced", "genes",
        ]
        writer = csv.DictWriter(buffer, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for entry in self:
            row = dict(entry)
            row["genes"] = ";".join(entry.get("genes", []))
            writer.writerow(row)
        text = buffer.getvalue()
        if path is not None:
            with open(path, "w", encoding="utf-8", newline="") as handle:
                handle.write(text)
        return text

    def to_text(self) -> str:
        """Return a compact, human-readable table."""
        if not self:
            return "No changes were made."
        lines = [
            f"{'Reaction':<24} {'Compartment':<12} {'Old bounds':<22} "
            f"{'New bounds':<22} Reason",
            "-" * 110,
        ]
        for entry in self:
            old = f"[{entry['old_lb']:g}, {entry['old_ub']:g}]"
            new = f"[{entry['new_lb']:g}, {entry['new_ub']:g}]"
            lines.append(
                f"{entry['reaction_id']:<24} {entry['compartment']:<12} "
                f"{old:<22} {new:<22} {entry['reason']}"
            )
        return "\n".join(lines)

    def __str__(self):  # pragma: no cover - convenience only
        return self.to_text()


# ----------------------------------------------------------------------
# model inspection helpers
# ----------------------------------------------------------------------
def _reaction_id(model, reaction_index: int) -> str:
    try:
        return str(model.rxns[reaction_index])
    except Exception:
        return str(reaction_index)


def _reaction_name(model, reaction_index: int) -> str:
    try:
        names = getattr(model, "rxnNames", None) or []
        if reaction_index < len(names):
            return str(names[reaction_index])
    except Exception:
        pass
    return ""


def _bound_pair(bounds: Sequence[float]):
    values = list(bounds)
    if len(values) != 2:
        raise ValueError("bounds must contain exactly two values")
    return float(values[0]), float(values[1])


def metabolite_compartment(model, metabolite_index: int) -> str:
    """Return the compartment abbreviation of one metabolite."""
    compartments = getattr(model, "metCompartments", None)
    if compartments is not None and metabolite_index < len(compartments):
        value = str(compartments[metabolite_index]).strip()
        if value:
            return value

    try:
        met_id = str(model.mets[metabolite_index])
    except Exception:
        return "?"

    if "_" not in met_id:
        return "?"

    suffix = met_id.rsplit("_", 1)[-1]
    known = list(getattr(model, "compartments", None) or [])
    # Only treat the suffix as a compartment if it is a known compartment or
    # no compartment list is available (the usual BiGG situation).
    if not known or suffix in known:
        return suffix
    return "?"


def reaction_compartment(model, reaction_index: int) -> str:
    """Return the compartment(s) touched by a reaction, e.g. ``"c"``."""
    try:
        matrix = model.S.toarray() if hasattr(model.S, "toarray") else np.asarray(model.S)
        column = matrix[:, int(reaction_index)]
    except Exception:
        return "?"
    compartments = []
    for met_index in np.flatnonzero(np.abs(column) > 1e-12):
        compartment = metabolite_compartment(model, int(met_index))
        if compartment and compartment not in compartments:
            compartments.append(compartment)
    return ",".join(compartments) if compartments else "?"


def compartment_names(compartment: str) -> str:
    """Translate an abbreviation such as ``"c"`` to a descriptive name."""
    if not compartment:
        return ""
    names = []
    for part in str(compartment).split(","):
        part = part.strip()
        names.append(COMPARTMENT_NAMES.get(part, part))
    return ", ".join(names)


def reaction_formula(model, reaction_index: int) -> str:
    """Build a reaction equation from ``S`` and ``metFormulas``.

    Metabolite formulae are used when available; otherwise metabolite IDs are
    used.  Exchange and demand reactions are rendered as ``substrate -->``.
    """
    try:
        matrix = model.S.toarray() if hasattr(model.S, "toarray") else np.asarray(model.S)
        column = matrix[:, int(reaction_index)]
    except Exception:
        return ""

    met_formulas = getattr(model, "metFormulas", None) or []
    reactants = []
    products = []
    for met_index, coefficient in enumerate(column):
        coefficient = float(coefficient)
        if abs(coefficient) <= 1e-12:
            continue
        formula = ""
        if met_index < len(met_formulas):
            formula = str(met_formulas[met_index]).strip()
        if not formula:
            try:
                formula = str(model.mets[met_index])
            except Exception:
                formula = f"met{met_index}"
        if abs(abs(coefficient) - 1.0) > 1e-9:
            formula = f"{abs(coefficient):g} {formula}"
        if coefficient < 0:
            reactants.append(formula)
        else:
            products.append(formula)

    if not reactants and not products:
        return ""
    return f"{' + '.join(reactants)} --> {' + '.join(products)}".strip()


def gpr_genes(rule) -> list:
    """Return the gene identifiers appearing in a GPR rule."""
    if rule is None:
        return []
    tokens = re.findall(r"[^\s()]+", str(rule))
    return [token for token in tokens if token not in ("and", "or")]


def write_change_log(log: ChangeLog, path: str):
    """Write ``log`` to JSON, CSV or plain text based on the suffix."""
    suffix = Path(path).suffix.lower()
    if suffix == ".json":
        log.to_json(path)
    elif suffix == ".csv":
        log.to_csv(path)
    else:
        Path(path).write_text(log.to_text(), encoding="utf-8")
    return path
