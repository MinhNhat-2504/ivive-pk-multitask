"""State which version of the deposit each paper used.

The concept DOI resolves to whatever version is current, which is what a
citation should do, but on its own it does not say what the authors actually
worked from. Version 2 is the one both papers were written against, so the
reference entries and the data statements name it alongside the DOI. If a
version 3 ever appears the citation still resolves and the paper still says
what it used.

Checked against the Zenodo API before writing: record 23471468 is version v2,
its concept DOI is 23462435, and all twelve files match the local bundle by
checksum.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPERS = ROOT / "paper"

EDITS = {
    "paper_A_supporting_information.md": [
        ("Archived at Zenodo, DOI: 10.5281/zenodo.23462435 (CC BY 4.0).",
         "Archived at Zenodo, DOI: 10.5281/zenodo.23462435, version v2 "
         "(CC BY 4.0). The DOI resolves to the current version; v2 is the one "
         "this paper was written against."),
    ],
    "paper_C_supporting_information.md": [
        ("(DOI 10.5281/zenodo.23462435).",
         "(DOI 10.5281/zenodo.23462435, version v2)."),
        ("DOI 10.5281/zenodo.23462435.",
         "DOI 10.5281/zenodo.23462435, version v2."),
    ],
    "paper_A_cover_letter.md": [
        ("are deposited at Zenodo under CC BY 4.0, DOI 10.5281/zenodo.23462435.",
         "are deposited at Zenodo under CC BY 4.0, DOI 10.5281/zenodo.23462435 "
         "(version v2)."),
    ],
    "paper_C_cover_letter.md": [
        ("DOI 10.5281/zenodo.23462435. Analysis code, including scripts that regenerate",
         "DOI 10.5281/zenodo.23462435, version v2. Analysis code, including scripts "
         "that regenerate"),
    ],
}

#: Paper C's deposit reference appears twice -- once in the rendered list and
#: once in the keyed block that `renumber_refs.py` regenerates it from -- so
#: both copies are replaced. Paper A carries the DOI in prose instead.
REFERENCE_EDITS = {
    "paper_C_target_scale.md": [
        ("Zenodo, 2026. DOI: 10.5281/zenodo.23462435.",
         "Zenodo, 2026, version v2. DOI: 10.5281/zenodo.23462435."),
    ],
    "paper_A_benchmark_defects.md": [
        ("are archived at Zenodo, DOI: 10.5281/zenodo.23462435 (CC BY 4.0)",
         "are archived at Zenodo, DOI: 10.5281/zenodo.23462435, version v2 "
         "(CC BY 4.0)"),
    ],
}


def apply(path: Path, edits: list[tuple[str, str]], every: bool = False) -> int:
    """Apply the edits. `every` replaces all copies rather than requiring one.

    A reference that `renumber_refs.py` manages exists twice: in the rendered
    list and in the keyed block it is regenerated from. Changing one and not
    the other would be silently undone on the next renumbering run.
    """
    if not path.exists():
        return 0
    text = path.read_text(encoding="utf-8")
    applied = 0
    for old, new in edits:
        if new in text and old not in text:
            continue
        count = text.count(old)
        if count == 0 or (not every and count != 1):
            continue
        text = text.replace(old, new)
        applied += count
    path.write_text(text, encoding="utf-8")
    return applied


def main() -> None:
    total = 0
    for name, edits in EDITS.items():
        count = apply(PAPERS / name, edits)
        print(f"  {name}: {count} statement(s) annotated")
        total += count
    for name, edits in REFERENCE_EDITS.items():
        count = apply(PAPERS / name, edits, every=True)
        print(f"  {name}: {count} copy/copies annotated")
        total += count
    print(f"\n{total} statement(s) now name the deposit version")


if __name__ == "__main__":
    main()
