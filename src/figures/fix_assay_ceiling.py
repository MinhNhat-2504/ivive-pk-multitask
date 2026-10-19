"""Retire a clip of my own invention from both manuscripts.

The fold error in unbound fraction was computed after clipping predictions at
99.99% bound, which is a number with nothing behind it. Where a model predicts
complete binding, fu is zero and the fold error is infinite, so some clip is
needed -- but the one to use is the assay's own declared ceiling of 99.95%,
beyond which no measurement could have reported anything.

The difference is not cosmetic for models that saturate. Clipping at 99.99
rather than at the declared ceiling turned a ridge model's median fold error
from 20 into 30 and a neural network's from 26.5 into 36.3, and the inflation
belonged to the clip. It also manufactured a gap between the two statistics for
the gradient-boosted model, 8.3 against 8.8, which disappears once the clip is
principled: both are 8.3.

So paper A loses the 8.8 it gained yesterday and keeps the 35-of-37 count,
which is real and is the stronger statement. Paper C's Table 2 is recomputed
and gains a column for how many compounds each model pushes to the ceiling,
since that is what made the figures clip-sensitive in the first place and is
itself the saturating end of the compression.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER_A = ROOT / "paper" / "paper_A_benchmark_defects.md"
PAPER_C = ROOT / "paper" / "paper_C_target_scale.md"

A_EDITS = [
    (
        """factor of 8.3 in 35 of the 37 compounds bound above 99% — 8.8-fold counting
both directions — because models trained""",
        """factor of 8.3 in 35 of the 37 compounds bound above 99%, because models trained""",
    ),
    (
        """is a median over-estimate of 8.3-fold, and a median fold error of 8.8, in the quantity that every well-stirred clearance calculation consumes.""",
        """is a median error of 8.3-fold in the quantity that every well-stirred clearance calculation consumes.""",
    ),
    (
        """37 cases, by a median factor of 8.3**; counting errors in both directions the
median fold error is 8.8. We give both because the first is the systematic
finding and the second is what a user of the number experiences. Because hepatic""",
        """37 cases, by a median factor of 8.3**. Counting errors in both directions gives
the same figure, because only two compounds fall the other way. Predictions are
clipped at the assay's declared ceiling of 99.95% bound, beyond which no
measurement could have reported anything; two of the 37 reach it. Because
hepatic""",
    ),
]

C_EDITS = [
    (
        """| Family | Spearman | Predicted range | Bias, bound < 50% | Bias, bound > 99% | Median $f_u$ over-estimate | Median $f_u$ fold error | Over-estimated |
|---|---|---|---|---|---|---|---|
| Ridge | 0.461 | 46.5 | +47.4 | −7.9 | 8.2 | **30.0** | 23 / 37 |
| SVR | 0.554 | 24.3 | +60.5 | −7.8 | 11.2 | 11.2 | 37 / 37 |
| $k$-NN | 0.305 | 73.1 | +47.0 | −15.5 | 17.9 | 17.9 | 37 / 37 |
| Random forest | 0.662 | 59.2 | +40.7 | −7.0 | 10.1 | 10.1 | 37 / 37 |
| Gradient boosting | 0.682 | 71.3 | +35.9 | −6.4 | 8.3 | **8.8** | 35 / 37 |
| Neural network | 0.466 | 68.8 | +28.1 | −19.4 | 26.5 | **36.3** | 31 / 37 |""",
        """| Family | Spearman | Predicted range | Bias, bound < 50% | Bias, bound > 99% | Median $f_u$ over-estimate | Median $f_u$ fold error | Over-estimated | At ceiling |
|---|---|---|---|---|---|---|---|---|
| Ridge | 0.461 | 46.5 | +47.4 | −7.9 | 8.2 | **20.0** | 23 / 37 | **14** |
| SVR | 0.554 | 24.3 | +60.5 | −7.8 | 11.2 | 11.2 | 37 / 37 | 0 |
| $k$-NN | 0.305 | 73.1 | +47.0 | −15.5 | 17.9 | 17.9 | 37 / 37 | 0 |
| Random forest | 0.662 | 59.2 | +40.7 | −7.0 | 10.1 | 10.1 | 37 / 37 | 0 |
| Gradient boosting | 0.682 | 71.3 | +35.9 | −6.4 | 8.3 | 8.3 | 35 / 37 | 2 |
| Neural network | 0.466 | 68.8 | +28.1 | −19.4 | 26.5 | 26.5 | 31 / 37 | **6**  |""",
    ),
    (
        """fraction: the median of $\\hat f_u / f_u$, which measures the systematic
over-estimate and lets opposite errors cancel, and the median of
$\\max(r, 1/r)$, which is what "wrong by a factor of" means. They coincide only
where every compound falls the same way, and the last column says how many do.""",
        """fraction: the median of $\\hat f_u / f_u$, which measures the systematic
over-estimate and lets opposite errors cancel, and the median of
$\\max(r, 1/r)$, which is what "wrong by a factor of" means. The last two columns
say how many compounds are over-estimated and how many the model pushes to or
past the assay's declared ceiling of 99.95% bound, at which $f_u$ is clipped —
a prediction of complete binding would otherwise give an infinite fold error,
and clipping anywhere tighter than the ceiling inflates the figure by an amount
the clip chooses rather than the model.""",
    ),
    (
        """The two $f_u$ columns part company exactly where a family's errors are not all
in one direction. Ridge over-estimates $f_u$ for 23 of the 37 compounds and
under-estimates it badly for the rest, which a signed median hides and a fold
error does not; reporting only the former would have understated its error by a
factor of nearly four. We note this because the figure of 8.3 that appears in
the companion paper¹³ is the signed statistic, whose fold-error
counterpart for the same model is 8.8.""",
        """The two $f_u$ columns part company for one family only. Ridge over-estimates
$f_u$ for 23 of the 37 compounds and places 14 of them at or beyond the assay's
ceiling, which a signed median hides and a fold error does not: 8.2 against
20.0. For the other five the two agree to the digits shown, so the figure of 8.3
reported for a gradient-boosted model in the companion paper[[prior]] stands
whichever statistic is meant.

The ceiling column is worth reading alongside the predicted range. A model that
compresses can do so by refusing the extremes, as SVR does with a range of 24.3
percentage points and no prediction at the ceiling, or by pushing highly bound
compounds past the point where the assay can distinguish them, as ridge does
with 14. Both destroy $f_u$; only the first is visible in a plot of percentage
points.""",
    ),
]


def apply(path: Path, edits: list[tuple[str, str]]) -> None:
    """Apply the edits, skipping any already present.

    The two files are patched in one run and the second failed once while the
    first had been written, which left the script unrunnable: paper A's
    anchors were gone. Skipping an edit whose result is already in the file
    makes a second run a no-op rather than an error.
    """
    text = path.read_text(encoding="utf-8")
    applied = 0
    for old, new in edits:
        if new in text and old not in text:
            continue
        if text.count(old) != 1:
            raise SystemExit(f"{path.name}: anchor not unique or missing:\n  {old[:90]}")
        text = text.replace(old, new, 1)
        applied += 1
    path.write_text(text, encoding="utf-8")
    print(f"{path.name}: {applied} passage(s) corrected, {len(edits) - applied} already done")


def main() -> None:
    apply(PAPER_A, A_EDITS)
    apply(PAPER_C, C_EDITS)


if __name__ == "__main__":
    main()
