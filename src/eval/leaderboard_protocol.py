"""Variance under the leaderboard's actual protocol, versus under resampling.

A correction. Earlier work in this project measured seed-to-seed variance with
`ADME(...).get_split(method="scaffold", seed=s)`, which draws a fresh split --
including a fresh *test* set -- on every seed. The leaderboard does not do
that. `admet_group` ships a **fixed** test set (243 molecules for
Clearance_Hepatocyte_AZ, identical hash on every seed) and the five seeds
repartition only train and validation, 848/122.

So the two protocols measure different things:

    leaderboard   fixed test; variance from the train partition and the learner
    resampled     fresh test each seed; adds the variance of which molecules
                  are being scored at all

Both are reported here against the same model and feature set, because the gap
between them is what says how much a published leaderboard number would move
on a different sample of molecules.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, PK_ENDPOINTS
from src.features.descriptors import featurise
from src.models.baseline import ENDPOINT_SPEC, LGBM_PARAMS, RESULTS_ROOT, _forward, _inverse

log = logging.getLogger(__name__)
N_DRAWS = 20_000
PROTOCOL_SEEDS = 5


def run_endpoint(group: admet_group, tdc_name: str, key: str, seeds: list[int]) -> pd.DataFrame:
    """Leaderboard protocol: fixed test set, train partition varies by seed."""
    spec = ENDPOINT_SPEC[key]
    benchmark = group.get(tdc_name)
    test = benchmark["test"]

    x_test, ok_test = featurise(test["Drug"].tolist())
    y_test = test["Y"].to_numpy()[ok_test]

    rows = []
    for seed in seeds:
        train, _valid = group.get_train_valid_split(
            benchmark=tdc_name, split_type="default", seed=seed
        )
        x_train, ok_train = featurise(train["Drug"].tolist())
        y_train = train["Y"].to_numpy()[ok_train]

        model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
        model.fit(x_train, _forward(y_train, spec["log"]))
        prediction = _inverse(model.predict(x_test), spec["log"])

        value = (
            float(spearmanr(y_test, prediction).statistic)
            if spec["metric"] == "spearman"
            else float(mean_absolute_error(y_test, prediction))
        )
        rows.append({"endpoint": key, "seed": seed, "metric": spec["metric"], "value": value})
        log.info("%s seed %2d: %.4f", key, seed, value)

    return pd.DataFrame(rows)


def null_gap(values: np.ndarray, rng: np.random.Generator) -> dict:
    """Apparent difference between two disjoint 5-seed runs of the same model."""
    gaps = np.empty(N_DRAWS)
    for i in range(N_DRAWS):
        picked = rng.choice(values, size=2 * PROTOCOL_SEEDS, replace=False)
        gaps[i] = picked[:PROTOCOL_SEEDS].mean() - picked[PROTOCOL_SEEDS:].mean()
    absolute = np.abs(gaps)
    return {
        "null_gap_median": float(np.median(absolute)),
        "null_gap_p95": float(np.percentile(absolute, 95)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 26)))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    group = admet_group(path=str(DATA_ROOT / "raw"))
    frames = [run_endpoint(group, name, key, args.seeds) for name, key in PK_ENDPOINTS.items()]
    per_seed = pd.concat(frames, ignore_index=True)

    rng = np.random.default_rng(0)
    rows = []
    resampled = pd.read_csv(RESULTS_ROOT / "baseline_lgbm_25_per_seed.csv")
    for (endpoint, metric), group_rows in per_seed.groupby(["endpoint", "metric"]):
        fixed = group_rows["value"].to_numpy()
        other = resampled.loc[resampled["endpoint"] == endpoint, "value"].to_numpy()
        rows.append(
            {
                "endpoint": endpoint,
                "metric": metric,
                "leaderboard_mean": fixed.mean(),
                "leaderboard_sd": fixed.std(ddof=1),
                **{f"leaderboard_{k}": v for k, v in null_gap(fixed, rng).items()},
                "resampled_mean": other.mean(),
                "resampled_sd": other.std(ddof=1),
                "sd_ratio": other.std(ddof=1) / fixed.std(ddof=1),
            }
        )
    summary = pd.DataFrame(rows)

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "leaderboard_protocol_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / "leaderboard_protocol_summary.csv", index=False)

    print("\n=== Leaderboard protocol (fixed test) vs resampled splits ===")
    print(summary.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
