"""Name the fu statistic correctly in both manuscripts.

Both papers quote "a median error of 8.3-fold" for the unbound fraction. That
number is the median of predicted/true, which is a signed ratio: it measures how
far fu is typically *over*-estimated and lets a compound in the other direction
cancel one in this direction. The fold error, max(r, 1/r), is what "wrong by a
factor of" means, and for this model it is 8.788 rather than 8.298.

The gap is small here because 35 of 37 compounds are over-estimated, which is
itself worth stating and is a stronger claim than the one the papers currently
make. It is not small everywhere: for a ridge model the signed median is 8.23
and the fold error 30.0, because only 23 of 37 fall the same way.

So the number stays and its name is corrected, the count is added, and the fold
error is given alongside where the claim is made in full.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
LETTER = ROOT / "paper" / "paper_A_cover_letter.md"

PAPER_EDITS = [
    # Abstract.
    ("""factor of 8.3 for the 29% of compounds bound above 99%, because models trained""",
     """factor of 8.3 in 35 of the 37 compounds bound above 99% — 8.8-fold counting
both directions — because models trained"""),
    # Contribution list.
    ("""binding value is wrong by a median factor of 8.3 for the compounds that matter""",
     """binding value over-states the unbound fraction by a median factor of 8.3 for
the compounds that matter"""),
    # Figure 3 caption.
    ("""a median error of 8.3-fold in the quantity that every well-stirred clearance calculation consumes.""",
     """a median over-estimate of 8.3-fold, and a median fold error of 8.8, in the quantity that every well-stirred clearance calculation consumes."""),
    # Section 4.5, where the claim is made in full.
    ("""predicted *f*u is 0.0617: **a median error of 8.3-fold**. Because hepatic""",
     """predicted *f*u is 0.0617. The unbound fraction is **over-estimated in 35 of the
37 cases, by a median factor of 8.3**; counting errors in both directions the
median fold error is 8.8. We give both because the first is the systematic
finding and the second is what a user of the number experiences. Because hepatic"""),
    # Recommendation.
    ("""expected to be highly bound. Section 4.5 puts the median error at 8.3-fold in""",
     """expected to be highly bound. Section 4.5 puts the median over-estimate at 8.3-fold in"""),
]

LETTER_EDITS = [
    ("""a median factor of 8.3 for the 29% of compounds bound above 99%, because a""",
     """a median factor of 8.3 in 35 of the 37 compounds bound above 99%, because a"""),
]


def apply(path: Path, edits: list[tuple[str, str]]) -> None:
    text = path.read_text(encoding="utf-8")
    for old, new in edits:
        if text.count(old) != 1:
            raise SystemExit(f"{path.name}: anchor not unique or missing:\n  {old[:80]}")
        text = text.replace(old, new, 1)
    path.write_text(text, encoding="utf-8")
    print(f"{path.name}: {len(edits)} passage(s) corrected")


def main() -> None:
    apply(PAPER, PAPER_EDITS)
    apply(LETTER, LETTER_EDITS)


if __name__ == "__main__":
    main()
