"""Point both papers at the deposit's concept DOI rather than at version one.

Zenodo issues two identifiers. The version DOI, 10.5281/zenodo.23462436, is
fixed to the deposit as it stood on 19 October 2026. The concept DOI,
10.5281/zenodo.23462435, resolves to whichever version is current.

Both papers were citing the version DOI. The deposit is about to gain the three
plasma protein binding species files, so a reader following the citation would
arrive at a version that does not contain the data the second paper rests on.
The concept DOI is the right one to cite for a dataset that will be versioned,
and the version is stated alongside it so the exact contents remain recoverable.

Checked against the Zenodo API before writing: `conceptdoi` is 23462435 and
`doi` is 23462436, so the two are not being confused.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPERS = ROOT / "paper"

VERSION_DOI = "10.5281/zenodo.23462436"
CONCEPT_DOI = "10.5281/zenodo.23462435"

FILES = [
    "paper_A_benchmark_defects.md",
    "paper_A_supporting_information.md",
    "paper_A_cover_letter.md",
    "paper_C_target_scale.md",
    "paper_C_supporting_information.md",
    "paper_C_cover_letter.md",
]


def main() -> None:
    total = 0
    for name in FILES:
        path = PAPERS / name
        if not path.exists():
            print(f"  {name}: absent, skipped")
            continue
        text = path.read_text(encoding="utf-8")
        count = text.count(VERSION_DOI)
        if not count:
            print(f"  {name}: nothing to change")
            continue
        path.write_text(text.replace(VERSION_DOI, CONCEPT_DOI), encoding="utf-8")
        print(f"  {name}: {count} occurrence(s) repointed")
        total += count

    tool = ROOT / "src" / "figures" / "set_doi.py"
    if tool.exists():
        text = tool.read_text(encoding="utf-8")
        if VERSION_DOI in text:
            tool.write_text(text.replace(VERSION_DOI, CONCEPT_DOI), encoding="utf-8")
            print("  set_doi.py: default repointed")

    print(f"\n{total} citation(s) now use the concept DOI {CONCEPT_DOI}")


if __name__ == "__main__":
    main()
