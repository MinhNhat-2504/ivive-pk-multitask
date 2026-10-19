"""A negative control for the target-scale mechanism, on a second endpoint.

The claim of `target_scale.py` is not that transformations help. It is specific:
a model should be fitted on the scale in which its consumer measures error, and
plasma protein binding is fitted wrongly because its label is a percentage while
its consumer, the well-stirred model, divides by fu and therefore cares about
relative error.

That claim is falsifiable, and this is the test. It predicts something about
endpoints where the convention already matches, not only about the one where it
does not. Intrinsic clearance is consumed by the same well-stirred model, and at
the extraction ratios this assay actually reaches -- median 0.022 on the 541
compounds with both a measured CLint and a measured binding value -- hepatic
clearance is near-linear in CLint, so the consumer's error in CL_h is
approximately the relative error in CLint. The appropriate scale is therefore
log CLint, which is already the field's convention and already what this project
uses.

So the predictions, recorded before running:

    1.  Fitting raw CLint instead of log CLint will be much worse in relative
        error, by a factor comparable to the PPB effect. The mechanism is the
        same -- squared error on a raw positive quantity prices a discrepancy by
        its absolute size, and the consumer prices it by its ratio.
    2.  Unlike PPB, no third scale will beat log CLint, because log already *is*
        the consumed scale. In particular the two scales that rescued PPB have
        no analogue here: there is no bounded complement to take.
    3.  Consequently the TDC leaderboard metric for this endpoint is not
        misaligned the way PPBR_AZ's is, and this should be visible as raw
        fitting losing on both the leaderboard metric and the relative error,
        rather than the two disagreeing.

If prediction 1 fails, the mechanism is wrong. If prediction 2 fails -- if some
further transform gains as much here as logit gained on PPB -- then the effect
is about transformation in general and not about matching the consumer, and the
paper's explanation is wrong even though its PPB result stands.

The benchmark metric for `Clearance_Microsome_AZ` is Spearman, so it is reported
alongside. The endpoint is a single species, unlike its hepatocyte counterpart.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr, wilcoxon
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)

SEEDS = list(range(1, 26))
#: The assay's declared working range, from the deposition's own description.
CLINT_FLOOR, CLINT_CEILING = 3.0, 150.0

#: Scales for a positive, unbounded-above quantity. `sqrt` and `reciprocal` are
#: included as foils for prediction 2: they are transformations, so if the
#: mechanism were "transformations help" they should help, and if the mechanism
#: is "match the consumer" they should not beat log.
SCALES = {
    "raw": (
        lambda y: np.asarray(y, dtype=float),
        lambda p: np.clip(p, 1e-3, None),
    ),
    "log10": (
        lambda y: np.log10(np.clip(y, 1e-6, None)),
        lambda p: 10.0 ** np.clip(p, -6.0, 6.0),
    ),
    "sqrt": (
        lambda y: np.sqrt(np.clip(y, 0.0, None)),
        lambda p: np.clip(p, 0.0, None) ** 2,
    ),
    "reciprocal": (
        lambda y: 1.0 / np.clip(y, 1e-3, None),
        lambda p: 1.0 / np.clip(p, 1e-6, None),
    ),
}


def score(truth: np.ndarray, pred: np.ndarray) -> dict:
    """Relative error, which is what the consumer of CLint sees, plus the
    leaderboard metric, which is a rank correlation."""
    fold = np.maximum(pred / truth, truth / pred)
    # The regime where the label is a detection-limit flag is excluded from the
    # fold-error summary: a value of 3.0 means "below 3", so a ratio against it
    # is not a measurement error. See paper A section 4.2.
    interior = (truth > CLINT_FLOOR) & (truth < CLINT_CEILING)
    return {
        "spearman": float(spearmanr(truth, pred).statistic),
        "mae_log10": float(np.mean(np.abs(np.log10(pred) - np.log10(truth)))),
        "fold_median": float(np.median(fold)),
        "fold_median_interior": float(np.median(fold[interior])),
        "frac_within_2fold": float(np.mean(fold <= 2.0)),
        "frac_within_2fold_interior": float(np.mean(fold[interior] <= 2.0)),
        "n_interior": int(interior.sum()),
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("Clearance_Microsome_AZ")

    x_train, ok_train = featurise(benchmark["train_val"]["Drug"].tolist())
    y_train = benchmark["train_val"]["Y"].to_numpy()[ok_train]
    x_test, ok_test = featurise(benchmark["test"]["Drug"].tolist())
    y_test = benchmark["test"]["Y"].to_numpy()[ok_test]

    interior = (y_test > CLINT_FLOOR) & (y_test < CLINT_CEILING)
    log.info(
        "train %d, test %d (%d interior to the declared range, %d on a bound)",
        len(y_train), len(y_test), int(interior.sum()), int((~interior).sum()),
    )

    rows = []
    for seed in SEEDS:
        for name, (forward, inverse) in SCALES.items():
            model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
            model.fit(x_train, forward(y_train))
            pred = inverse(model.predict(x_test))
            rows.append({"seed": seed, "scale": name, **score(y_test, pred)})

    per_seed = pd.DataFrame(rows)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "second_endpoint_per_seed.csv", index=False)

    metrics = [c for c in per_seed.columns if c not in {"seed", "scale", "n_interior"}]
    summary = per_seed.groupby("scale")[metrics].agg(["mean", "std"])

    print("\n=== Clearance_Microsome_AZ, target scale, 25 seeds ===")
    for metric in metrics:
        print(f"\n{metric}")
        for name in SCALES:
            cell = summary.loc[name, metric]
            print(f"  {name:>11s}  {cell['mean']:9.4f} +/- {cell['std']:.4f}")

    print("\n=== Paired against `log10`, the consumed scale ===")
    wide = per_seed.pivot(index="seed", columns="scale", values=metrics)
    for metric in metrics:
        print(f"\n{metric}")
        for name in ("raw", "sqrt", "reciprocal"):
            delta = (wide[(metric, name)] - wide[(metric, "log10")]).to_numpy()
            p = wilcoxon(delta).pvalue if np.ptp(delta) > 0 else float("nan")
            print(
                f"  {name:>11s} - log10: {delta.mean():+9.4f} "
                f"[{np.percentile(delta, 2.5):+.4f}, {np.percentile(delta, 97.5):+.4f}]  "
                f"p={p:.2g}"
            )

    summary.to_csv(RESULTS_ROOT / "second_endpoint_summary.csv")
    print("\nA positive delta on a fold-error row means that scale is worse than log10.")


if __name__ == "__main__":
    main()
