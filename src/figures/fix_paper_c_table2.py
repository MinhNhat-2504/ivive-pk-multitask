"""Correct paper C's Table 2 to a genuine fold error, and say which is which.

Table 2 was headed "Median fu fold error" and carried the median of
predicted/true, a signed ratio in which a compound wrong in one direction
cancels one wrong in the other. For most families every compound falls the same
way and the two agree, but not for all: ridge over-estimates fu in only 23 of 37
cases and its fold error is 30.0 where the signed median is 8.2, and the neural
network's are 36.3 and 26.5. The table understated two of six rows.

Both statistics are now given, because they answer different questions. The
signed median is the systematic finding -- fu is over-estimated -- and the fold
error is what a user of the number experiences.

The abstract already quotes the paired 25-seed figures, which were computed as
fold errors and need no change. Section 3.5 now states why its percent column
differs from Table 2's: targets are standardised there so that scale-dependent
hyperparameters mean the same thing on all three scales, which changes the
percent fit slightly.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_C_target_scale.md"

EDITS = [
    (
        """**Table 2. Compression on the external set, percent-bound target.** Observed
range is 93.9 percentage points.

| Family | Spearman | Predicted range | Bias, bound < 50% | Bias, bound > 99% | Median $f_u$ fold error |
|---|---|---|---|---|---|
| Ridge | 0.461 | 46.5 | +47.4 | −7.9 | 8.2 |
| SVR | 0.554 | 24.3 | +60.5 | −7.8 | 11.2 |
| $k$-NN | 0.305 | 73.1 | +47.0 | −15.5 | 17.9 |
| Random forest | 0.662 | 59.2 | +40.7 | −7.0 | 10.1 |
| Gradient boosting | 0.682 | 71.3 | +35.9 | −6.4 | 8.3 |
| Neural network | 0.466 | 68.8 | +28.1 | −19.4 | 26.5 |""",
        """**Table 2. Compression on the external set, percent-bound target.** Observed
range is 93.9 percentage points. Two statistics are given for the unbound
fraction: the median of $\\hat f_u / f_u$, which measures the systematic
over-estimate and lets opposite errors cancel, and the median of
$\\max(r, 1/r)$, which is what "wrong by a factor of" means. They coincide only
where every compound falls the same way, and the last column says how many do.

| Family | Spearman | Predicted range | Bias, bound < 50% | Bias, bound > 99% | Median $f_u$ over-estimate | Median $f_u$ fold error | Over-estimated |
|---|---|---|---|---|---|---|---|
| Ridge | 0.461 | 46.5 | +47.4 | −7.9 | 8.2 | **30.0** | 23 / 37 |
| SVR | 0.554 | 24.3 | +60.5 | −7.8 | 11.2 | 11.2 | 37 / 37 |
| $k$-NN | 0.305 | 73.1 | +47.0 | −15.5 | 17.9 | 17.9 | 37 / 37 |
| Random forest | 0.662 | 59.2 | +40.7 | −7.0 | 10.1 | 10.1 | 37 / 37 |
| Gradient boosting | 0.682 | 71.3 | +35.9 | −6.4 | 8.3 | **8.8** | 35 / 37 |
| Neural network | 0.466 | 68.8 | +28.1 | −19.4 | 26.5 | **36.3** | 31 / 37 |""",
    ),
    (
        """narrower than the observed range throughout. A model's ranking ability does not
protect its $f_u$: gradient boosting ranks best and has the smallest $f_u$ error,
the neural network ranks poorly and has the largest, but $k$-NN ranks worst of
all and beats the neural network. The usable statement is the weaker one — **a
good position on a percent-bound leaderboard does not imply a usable unbound
fraction.**""",
        """narrower than the observed range throughout. A model's ranking ability does not
protect its $f_u$: gradient boosting ranks best and has the smallest fold error,
the neural network ranks poorly and has the largest, but $k$-NN ranks worst of
all and beats the neural network, while ridge ranks in the middle and has the
second worst. The usable statement is therefore the weaker one — **a good
position on a percent-bound leaderboard does not imply a usable unbound
fraction.**

The two $f_u$ columns part company exactly where a family's errors are not all
in one direction. Ridge over-estimates $f_u$ for 23 of the 37 compounds and
under-estimates it badly for the rest, which a signed median hides and a fold
error does not; reporting only the former would have understated its error by a
factor of nearly four. We note this because the figure of 8.3 that appears in
the companion paper[[prior]] is the signed statistic, whose fold-error
counterpart for the same model is 8.8.""",
    ),
    (
        """| Family | `percent` | `log_fu` | `logit` | Improvement | Test MAE, `percent` → best |""",
        """Targets are standardised in every cell of this comparison, so the `percent`
column differs slightly from Table 2, which uses the models as an unaware
practitioner would fit them.

| Family | `percent` | `log_fu` | `logit` | Improvement | Test MAE, `percent` → best |""",
    ),
]


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    for old, new in EDITS:
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique or missing:\n  {old[:90]}")
        text = text.replace(old, new, 1)
    PAPER.write_text(text, encoding="utf-8")
    print(f"paper C: {len(EDITS)} passage(s) corrected")


if __name__ == "__main__":
    main()
