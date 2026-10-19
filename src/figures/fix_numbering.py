"""Realign citation markers after a reference was removed, and relocate one.

Dropping the Butina entry shifted references 17-19 down to 16-18, but the
markers had already been written against the old numbering. This corrects them
and resolves a second problem the check surfaced: Bemis and Murcko is cited
only from the Supporting Information, so it does not belong in the main
reference list. Scaffold overlap is tabulated in S1.5 and discussed nowhere in
the main text, so the entry moves to an SI reference list rather than being
cited gratuitously in the paper to justify its own presence.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"

SUPERSCRIPT = {str(d): c for d, c in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹")}
ENTRY = re.compile(r"^\((\d+)\) ")


def to_superscript(number: int) -> str:
    return "".join(SUPERSCRIPT[d] for d in str(number))


def main() -> None:
    paper = PAPER.read_text(encoding="utf-8")

    # Markers written against the pre-removal numbering.
    for old, new in ((19, 18), (18, 17)):
        paper = paper.replace(to_superscript(old), to_superscript(new))

    # Move Bemis & Murcko out of the main list; it is cited only from the SI.
    lines = paper.split("\n")
    start = next(i for i, line in enumerate(lines) if line.startswith("(16) Bemis"))
    end = start + 1
    while end < len(lines) and lines[end].strip() and not ENTRY.match(lines[end]):
        end += 1
    moved = " ".join(line.strip() for line in lines[start:end])
    while end < len(lines) and not lines[end].strip():
        end += 1

    kept = []
    for line in lines[:start] + lines[end:]:
        match = ENTRY.match(line)
        if match and int(match.group(1)) > 16:
            number = int(match.group(1))
            line = line.replace(f"({number})", f"({number - 1})", 1)
        kept.append(line)
    PAPER.write_text("\n".join(kept), encoding="utf-8")

    # The SI now needs its own single-entry reference list.
    si = SI.read_text(encoding="utf-8")
    si = si.replace(to_superscript(17), to_superscript(1))
    si = si.rstrip() + (
        "\n\n## S6. References cited only in this Supporting Information\n\n"
        "(S1) " + moved.split(") ", 1)[1] + "\n\n"
        "All other references are numbered as in the main text.\n"
    )
    si = si.replace("framework¹ is ordinary", "framework^S1^ is ordinary")
    SI.write_text(si, encoding="utf-8")

    print("markers realigned to 17/18; Bemis & Murcko moved to the SI as S1")


if __name__ == "__main__":
    main()
