"""Fill the author block in both manuscript files.

Also corrects the status line, which said "Single author" throughout because
the draft began that way.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"

AFFILIATION = (
    "Faculty of Information Technology, Van Lang School of Technology, "
    "Van Lang University, 69/68 Dang Thuy Tram Street, Binh Loi Trung Ward, "
    "Ho Chi Minh City, Vietnam"
)

BLOCK = f"""**[FIRST AUTHOR FULL NAME]**ᵃ and **Nguyen Thi My Linh**ᵃ'\\*

ᵃ {AFFILIATION}

ORCID: [first author] 0009-0006-9477-9642; N.T.M. Linh 0000-0003-3990-312X

\\* Corresponding author: linh.ntm@vlu.edu.vn

`[TODO: first author's full name as it should appear, and preferred romanisation]`"""

SI_BLOCK = f"""**[FIRST AUTHOR FULL NAME]** and **Nguyen Thi My Linh**\\*

{AFFILIATION}

\\* Corresponding author: linh.ntm@vlu.edu.vn"""


def main() -> None:
    paper = PAPER.read_text(encoding="utf-8")
    paper = paper.replace("`[TODO: author block, affiliation, ORCID]`", BLOCK)
    # The draft described itself as single-author from the start.
    paper = paper.replace(
        "Target: *Journal of Chemical Information and Modeling*. Single\nauthor.",
        "Target: *Journal of Chemical Information and Modeling*.",
    )
    paper = paper.replace(
        "Target: *Journal of Chemical Information and Modeling*. Single author.",
        "Target: *Journal of Chemical Information and Modeling*.",
    )
    PAPER.write_text(paper, encoding="utf-8")

    si = SI.read_text(encoding="utf-8")
    si = si.replace("`[TODO: author block]`", SI_BLOCK)
    SI.write_text(si, encoding="utf-8")
    print("author block written to both files")


if __name__ == "__main__":
    main()
