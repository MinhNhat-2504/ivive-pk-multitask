"""Write the author names into both manuscript files.

Names are used exactly as the authors gave them. For a Vietnamese name an
English-language journal needs to be told which element is the surname for
indexing, so the submission note in the cover letter states it rather than
leaving a copy editor to guess.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"

FIRST = "Nhat Trinh Ngoc Minh"
SECOND = "Linh Nguyen Thi My"

EDITS = [
    ("**[FIRST AUTHOR FULL NAME]**ᵃ and **Nguyen Thi My Linh**",
     f"**{FIRST}**ᵃ and **{SECOND}**"),
    ("**[FIRST AUTHOR FULL NAME]** and **Nguyen Thi My Linh**",
     f"**{FIRST}** and **{SECOND}**"),
    ("ORCID: [first author] 0009-0006-9477-9642; N.T.M. Linh 0000-0003-3990-312X",
     f"ORCID: {FIRST} 0009-0006-9477-9642; {SECOND} 0000-0003-3990-312X"),
    ("`[TODO: first author's full name as it should appear, and preferred romanisation]`",
     ""),
]


def main() -> None:
    for path in (PAPER, SI):
        text = path.read_text(encoding="utf-8")
        applied = 0
        for old, new in EDITS:
            if old in text:
                text = text.replace(old, new)
                applied += 1
        path.write_text(text, encoding="utf-8")
        print(f"{path.name}: {applied} edit(s)")


if __name__ == "__main__":
    main()
