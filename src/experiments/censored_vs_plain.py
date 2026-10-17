"""Does modelling the detection limits help? Paired per split.

Both arms see the identical train/test partition on every seed, so the
difference is taken within a split and the split variance -- which Section 4
showed to be the dominant term -- cancels. An unpaired comparison at this
effect size would be uninterpretable.

Reported per endpoint:

* the paired mean difference with a bootstrap interval and a Wilcoxon signed-
  rank test over seeds;
* the same on the uncensored test rows alone. The censored test labels are
  ties, and a rank metric computed over them rewards agreeing with a flag
  rather than predicting a quantity, so the uncensored subset is the cleaner
  read even though the full test set is what the leaderboard scores.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import wilcoxon
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT, PK_ENDPOINTS
from src.eval.censoring import ASSAY_LIMITS
from src.features.descriptors import featurise
from src.models.baseline import ENDPOINT_SPEC, LGBM_PARAMS, RESULTS_ROOT, _forward, score
from src.models.censored import CensoredLGBM

log = logging.getLogger(__name__)
N_BOOT = 10_000
RNG_SEED = 0


def run_seed(dataset: ADME, key: str, seed: int) -> list[dict]:
    """Fit both arms on one split and score them on the same test rows."""
    spec = ENDPOINT_SPEC[key]
    lo, hi = ASSAY_LIMITS[key]
    split = dataset.get_split(method="scaffold", seed=seed)
    train = pd.concat([split["train"], split["valid"]], ignore_index=True)
    test = split["test"]

    x_train, ok_train = featurise(train["Drug"].tolist())
    x_test, ok_test = featurise(test["Drug"].tolist())
    y_train_raw = train["Y"].to_numpy()[ok_train]
    y_test_raw = test["Y"].to_numpy()[ok_test]

    y_train = _forward(y_train_raw, spec["log"])
    lo_s, hi_s = _forward(np.array([lo, hi]), spec["log"])

    plain = LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_train, y_train)
    tobit = CensoredLGBM(lo_s, hi_s, random_state=seed, **LGBM_PARAMS).fit(x_train, y_train)

    uncensored_test = ~(np.isclose(y_test_raw, lo) | np.isclose(y_test_raw, hi))
    rows = []
    for subset, mask in (("all test rows", np.ones_like(uncensored_test)), ("uncensored test rows", uncensored_test)):
        if mask.sum() < 10:
            continue
        rows.append(
            {
                "endpoint": key,
                "seed": seed,
                "subset": subset,
                "n_test": int(mask.sum()),
                "plain": score(y_test_raw[mask], plain.predict(x_test)[mask], spec["metric"]),
                "tobit": score(y_test_raw[mask], tobit.predict(x_test)[mask], spec["metric"]),
            }
        )
    return rows


def paired_summary(frame: pd.DataFrame, higher_is_better: bool) -> dict:
    """Paired difference over seeds, with a bootstrap interval and a signed-rank test."""
    diff = (frame["tobit"] - frame["plain"]).to_numpy()
    if not higher_is_better:
        diff = -diff  # express as improvement regardless of metric direction

    rng = np.random.default_rng(RNG_SEED)
    boot = np.array(
        [rng.choice(diff, size=diff.size, replace=True).mean() for _ in range(N_BOOT)]
    )
    try:
        p_value = float(wilcoxon(diff).pvalue)
    except ValueError:  # all differences identical
        p_value = float("nan")

    return {
        "n_seeds": int(diff.size),
        "mean_improvement": float(diff.mean()),
        "ci_lo": float(np.percentile(boot, 2.5)),
        "ci_hi": float(np.percentile(boot, 97.5)),
        "wins": int((diff > 0).sum()),
        "wilcoxon_p": p_value,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints", nargs="*", default=list(ASSAY_LIMITS))
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 26)))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    name_of = {key: name for name, key in PK_ENDPOINTS.items()}

    records = []
    for key in args.endpoints:
        dataset = ADME(name=name_of[key], path=str(DATA_ROOT / "raw"))
        for seed in args.seeds:
            rows = run_seed(dataset, key, seed)
            records.extend(rows)
            shown = {r["subset"]: (r["plain"], r["tobit"]) for r in rows}
            log.info("%s seed %2d: %s", key, seed, shown)

    per_seed = pd.DataFrame(records)

    summaries = []
    for (endpoint, subset), group in per_seed.groupby(["endpoint", "subset"]):
        higher = ENDPOINT_SPEC[endpoint]["metric"] == "spearman"
        summaries.append(
            {
                "endpoint": endpoint,
                "subset": subset,
                "metric": ENDPOINT_SPEC[endpoint]["metric"],
                "plain_mean": group["plain"].mean(),
                "tobit_mean": group["tobit"].mean(),
                **paired_summary(group, higher),
            }
        )
    summary = pd.DataFrame(summaries)

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "censored_vs_plain_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / "censored_vs_plain_summary.csv", index=False)

    print("\n=== Paired comparison: Tobit objective vs ordinary least squares ===")
    print(summary.round(4).to_string(index=False))
    print(
        "\nmean_improvement is signed so that positive always means Tobit is better."
        f"\nwins counts seeds where Tobit won, out of n_seeds."
    )


if __name__ == "__main__":
    main()
