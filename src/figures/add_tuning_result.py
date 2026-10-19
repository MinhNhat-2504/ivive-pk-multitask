"""Add the tuning comparison to paper C, and retire the limitation it answers.

The manuscript's models are untuned, which invites the objection that the
percent-bound target is simply under-served by hyperparameters chosen without it
in mind. The limitations section conceded the point. It no longer has to: both
scales were given the same sixteen-configuration search over the same grid,
both selected on validation mean absolute error, which is the benchmark's own
metric and what a leaderboard entrant would tune for.

The result goes in Section 3.4, where the leaderboard is discussed, rather than
in a new subsection, so nothing downstream has to be renumbered.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_C_target_scale.md"

ADDITION = """

**Tuning does not substitute for the scale.** A reviewer may reasonably suspect
that the percent-bound target is simply under-served by hyperparameters chosen
without it in mind. Both were therefore given the same search: sixteen
configurations over the same grid, including the untuned one, selected on
validation mean absolute error — the benchmark's metric, and what an entrant
would tune for. The search moved away from the default configuration for every
scale and every seed, and improved validation error in all of them.

**Table 5. Sixteen configurations per scale, selected on validation MAE, five
seeds under the official protocol.**

| Target | Test MAE, untuned → tuned | Median $f_u$ fold error, untuned → tuned |
|---|---|---|
| `percent` | 8.009 → 7.991 | 8.97 → 8.70 |
| `log_fu` | 7.698 → 7.611 | 2.45 → 2.65 |
| `logit` | 7.553 → **7.410** | 2.59 → 2.63 |

Tuning improves mean absolute error on every scale, and improves it most on the
transformed ones. It does not close the $f_u$ gap at all. A tuned percent model
reaches 7.991 and an $f_u$ error of 8.70; an **untuned** logit model reaches
7.553 and 2.59, which is better on both. Across the five seeds and both tuning
states the two sets of $f_u$ errors do not overlap: 6.98 to 10.34 for percent,
2.04 to 3.25 for the transformed scales.

This is what the argument of Section 3.3 predicts. The selection criterion is
mean absolute error on percent bound, which is nearly blind to the region where
$f_u$ varies, so searching under it cannot recover what the target scale
discards. The effect is visible in the other direction too: tuning on that
criterion slightly *degrades* the $f_u$ accuracy of a well-scaled model, taking
`log_fu` from 2.45-fold to 2.65 and its fraction within two-fold from 0.384 to
0.330."""

EDITS = [
    # The three tables after the insertion point shift by one. Renamed from the
    # bottom up: renaming 6 to 7 before 7 to 8 would make the second rename
    # match the table the first had just created. Sequential renumbering has
    # already sent two citations to the same marker once in this project.
    (
        """**Table 7. Median $f_u$ fold error on the external set, compounds bound above
99%, 25 seeds.**""",
        """**Table 8. Median $f_u$ fold error on the external set, compounds bound above
99%, 25 seeds.**""",
    ),
    (
        """**Table 6. Median $f_u$ fold error for compounds bound above 99%, and the
benchmark metric alongside.**""",
        """**Table 7. Median $f_u$ fold error for compounds bound above 99%, and the
benchmark metric alongside.**""",
    ),
    (
        """**Table 5. Target scale, gradient boosting, 25 seeds, trained on train-validation.**""",
        """**Table 6. Target scale, gradient boosting, 25 seeds, trained on train-validation.**""",
    ),
    (
        """Trained instead on the whole train-validation partition, which is not the
leaderboard's protocol but is the fairer comparison between models, the same""",
        ADDITION.strip() + """

Trained instead on the whole train-validation partition, which is not the
leaderboard's protocol but is the fairer comparison between models, the same""",
    ),
    # The limitation is answered.
    (
        """Our leaderboard positions are computed under the benchmark's protocol with an
untuned model and no validation-based early stopping, so they are a floor on what
the transformation is worth rather than an attempt at a competitive entry. We
have not re-run published entries on a transformed target, which would be the
direct test of whether the ranking changes throughout.""",
        """Our leaderboard positions are computed under the benchmark's protocol with a
small search and no validation-based early stopping, so they are a floor on what
the transformation is worth rather than an attempt at a competitive entry. We
have not re-run published entries on a transformed target, which would be the
direct test of whether the ranking changes throughout; Section 3.4 shows only
that a search on one scale does not reach another scale's untuned baseline.""",
    ),
]


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    applied = 0
    for old, new in EDITS:
        if new in text and old not in text:
            continue
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique or missing:\n  {old[:90]}")
        text = text.replace(old, new, 1)
        applied += 1
    PAPER.write_text(text, encoding="utf-8")
    print(f"paper C: {applied} passage(s) updated")


if __name__ == "__main__":
    main()
