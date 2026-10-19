"""Write the Zenodo DOI into both manuscript files.

The version DOI is used, not the concept DOI. The concept DOI resolves to
whichever version is newest, which is the wrong citation for a paper whose
argument is that provenance must be fixed and inspectable.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"

DOI = "10.5281/zenodo.23462435"

PAPER_TEXT = (
    "The verified external set, the hand-review files recording every accept, "
    "correct, reject and withhold judgement, and the AstraZeneca assays as "
    "deposited with species and measurement qualifier intact, are archived at "
    f"Zenodo, DOI: {DOI} (CC BY 4.0). The version DOI is cited rather than the "
    "concept DOI so that the citation cannot come to mean a different deposit."
)

SI_TEXT = f"Archived at Zenodo, DOI: {DOI} (CC BY 4.0)."


def main() -> None:
    paper = PAPER.read_text(encoding="utf-8")
    old = "`[TODO: Zenodo DOI for the dataset]`"
    if old not in paper:
        raise SystemExit("placeholder not found in the manuscript")
    PAPER.write_text(paper.replace(old, PAPER_TEXT), encoding="utf-8")

    si = SI.read_text(encoding="utf-8")
    old_si = "`[TODO: Zenodo DOI once deposited]`"
    if old_si not in si:
        raise SystemExit("placeholder not found in the SI")
    SI.write_text(si.replace(old_si, SI_TEXT), encoding="utf-8")

    print(f"DOI {DOI} written to both files")


if __name__ == "__main__":
    main()
