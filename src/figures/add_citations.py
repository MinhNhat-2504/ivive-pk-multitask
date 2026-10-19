"""Attach the citation markers the consistency check found missing.

`check_paper.py` reported nine references listed but never cited, and three
`[ref]` placeholders survived the split from the combined draft. Eight
references have a place in the text and are given one here; the Butina
clustering paper has none, because the cluster splits it supports belong to the
companion method paper, so it is removed rather than cited gratuitously.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"

PAPER_EDITS: list[tuple[str, str]] = [
    # RDKit (19), where canonical SMILES is first used.
    (
        "Records are matched by RDKit canonical SMILES",
        "Records are matched by RDKit¹⁹ canonical SMILES",
    ),
    # Temesi (7), the assay platform behind the deposition.
    (
        "deposition's own assay description, and separately",
        "deposition's own assay description — the assays are the high-throughput "
        "metabolic-stability platform described by Temesi et al.⁷ — and separately",
    ),
    # LightGBM (18), the one model the paper fits.
    (
        "Fitting one untuned model under both",
        "Fitting one untuned gradient-boosted model¹⁸ under both",
    ),
    # Tobin (15), at the first mention of the Tobit likelihood.
    (
        "that a Tobit likelihood modelling the censored labels changes nothing",
        "that a Tobit likelihood¹⁵ modelling the censored labels changes nothing",
    ),
    # Three [ref] placeholders left from the split.
    (
        "The source deposition is ChEMBL CHEMBL3301361 [ref], an AstraZeneca direct",
        "The source deposition is ChEMBL CHEMBL3301361,⁶ an AstraZeneca direct",
    ),
    (
        "and cites Di et al. [ref] on human microsomes",
        "and cites Di et al.⁸ on human microsomes",
    ),
    # Bowman & Benet (12), Wilson (13), Obach (9) and Gidal (14), where the
    # paper states that no repeatability estimate survives for these endpoints.
    (
        "The human assays carry one record per molecule, and `Half_Life_Obach` has no\n"
        "replicate structure by construction, being a compilation of one intravenous\n"
        "study per drug [ref].",
        "The human assays carry one record per molecule. Literature estimates exist "
        "but do not substitute: interlaboratory coefficients of variation for human "
        "hepatocyte intrinsic clearance span 8.5% to 99.8% across compounds compiled "
        "from fourteen sources,¹² and the scaling factors that would relate "
        "them vary several-fold between donors,¹³ so a borrowed ceiling "
        "carries no measurement of *this* assay. `Half_Life_Obach` has no replicate "
        "structure by construction, being a compilation of one intravenous study per "
        "drug,⁹ and terminal half-life is in any case design-dependent: the same "
        "topiramate data yield 28.8 or 82.8 hours depending only on how long sampling "
        "continued.¹⁴",
    ),
]

SI_EDITS: list[tuple[str, str]] = [
    # Bemis & Murcko (17), where the SI tabulates scaffold overlap.
    (
        "sharing a Bemis–Murcko framework is ordinary chemical space",
        "sharing a Bemis–Murcko framework¹⁷ is ordinary chemical space",
    ),
]

#: Supports nothing in this paper; the cluster splits it underpins are in the
#: companion method paper.
REMOVE = "(16) Butina, D."


def apply(path: Path, edits: list[tuple[str, str]]) -> None:
    text = path.read_text(encoding="utf-8")
    for old, new in edits:
        if old not in text:
            raise SystemExit(f"anchor not found in {path.name}:\n  {old[:90]}")
        text = text.replace(old, new, 1)
    path.write_text(text, encoding="utf-8")
    print(f"{path.name}: {len(edits)} edit(s) applied")


#: A reference entry starts a line as "(12) ". Prose inside an entry can also
#: open with a bracket -- the pinned RDKit entry does -- so the pattern has to
#: require a number and a space, not merely a bracket.
ENTRY = re.compile(r"^\((\d+)\) ")


def drop_reference(path: Path, marker: str) -> None:
    """Remove one reference entry and renumber those after it."""
    lines = path.read_text(encoding="utf-8").split("\n")
    start = next(i for i, line in enumerate(lines) if line.startswith(marker))
    end = start + 1
    while end < len(lines) and lines[end].strip() and not ENTRY.match(lines[end]):
        end += 1
    while end < len(lines) and not lines[end].strip():
        end += 1

    removed = int(ENTRY.match(lines[start]).group(1))
    renumbered = []
    for line in lines[:start] + lines[end:]:
        match = ENTRY.match(line)
        if match and int(match.group(1)) > removed:
            number = int(match.group(1))
            line = line.replace(f"({number})", f"({number - 1})", 1)
        renumbered.append(line)
    path.write_text("\n".join(renumbered), encoding="utf-8")
    print(f"{path.name}: removed reference {removed}, renumbered the rest")


def main() -> None:
    apply(PAPER, PAPER_EDITS)
    apply(SI, SI_EDITS)
    drop_reference(PAPER, REMOVE)


if __name__ == "__main__":
    main()
