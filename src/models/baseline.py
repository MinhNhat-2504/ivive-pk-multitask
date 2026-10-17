"""Single-task LightGBM baselines on the TDC scaffold splits.

Purpose is reproduction, not improvement: before any constrained model is worth
reporting, we need our own numbers for the five PK endpoints under the same
protocol the leaderboard uses, with the seed-to-seed spread made visible. That
spread is what decides whether a later gain is real, and the published
leaderboard reports it only sparsely.

Protocol follows the TDC ADMET leaderboard: scaffold split, seeds 1-5, metric
Spearman for the four regression-on-rank endpoints and MAE for PPBR.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT, PK_ENDPOINTS
from src.features.descriptors import featurise

log = logging.getLogger(__name__)

SEEDS = [1, 2, 3, 4, 5]
RESULTS_ROOT = Path(__file__).resolve().parents[2] / "results"

#: Metric used by the TDC leaderboard for each endpoint, and whether the target
#: is log-transformed before fitting. Clearance, half-life and volume of
#: distribution span orders of magnitude and are modelled on a log scale as is
#: standard in PK; PPBR is a bounded percentage and is modelled as measured.
ENDPOINT_SPEC: dict[str, dict] = {
    "cl_hep": {"metric": "spearman", "log": True},
    "cl_mic": {"metric": "spearman", "log": True},
    "ppbr": {"metric": "mae", "log": False},
    "vdss": {"metric": "spearman", "log": True},
    "half_life": {"metric": "spearman", "log": True},
}

LGBM_PARAMS = dict(
    n_estimators=600,
    learning_rate=0.03,
    num_leaves=31,
    min_child_samples=10,
    subsample=0.8,
    subsample_freq=1,
    colsample_bytree=0.6,
    reg_lambda=1.0,
    n_jobs=-1,
    verbose=-1,
)


def _forward(y: np.ndarray, use_log: bool) -> np.ndarray:
    """Map labels to the space the model is fitted in."""
    if not use_log:
        return y
    # Clearance and half-life are non-negative and occasionally exactly zero;
    # log1p keeps those rows instead of discarding a real measurement.
    return np.log1p(np.clip(y, 0.0, None))


def _inverse(y: np.ndarray, use_log: bool) -> np.ndarray:
    """Map predictions back to the measured scale."""
    return np.expm1(y) if use_log else y


def score(y_true: np.ndarray, y_pred: np.ndarray, metric: str) -> float:
    if metric == "spearman":
        return float(spearmanr(y_true, y_pred).statistic)
    if metric == "mae":
        return float(mean_absolute_error(y_true, y_pred))
    raise ValueError(f"unknown metric: {metric}")


def run_endpoint(tdc_name: str, key: str, seeds: list[int] = SEEDS) -> pd.DataFrame:
    """Fit and evaluate one endpoint across seeds. One row per seed."""
    spec = ENDPOINT_SPEC[key]
    dataset = ADME(name=tdc_name, path=str(DATA_ROOT / "raw"))
    rows = []

    for seed in seeds:
        split = dataset.get_split(method="scaffold", seed=seed)
        # Validation rows are folded into training: no early stopping or tuning
        # happens here, so holding them out would only shrink the fit for no
        # methodological gain. The test fold stays untouched.
        train = pd.concat([split["train"], split["valid"]], ignore_index=True)
        test = split["test"]

        x_train, ok_train = featurise(train["Drug"].tolist())
        x_test, ok_test = featurise(test["Drug"].tolist())
        y_train = train["Y"].to_numpy()[ok_train]
        y_test = test["Y"].to_numpy()[ok_test]

        model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
        model.fit(x_train, _forward(y_train, spec["log"]))
        y_pred = _inverse(model.predict(x_test), spec["log"])

        value = score(y_test, y_pred, spec["metric"])
        rows.append(
            {
                "endpoint": key,
                "seed": seed,
                "metric": spec["metric"],
                "value": value,
                "n_train": len(y_train),
                "n_test": len(y_test),
            }
        )
        log.info("%s seed %d: %s = %.4f", key, seed, spec["metric"], value)

    return pd.DataFrame(rows)


def summarise(per_seed: pd.DataFrame) -> pd.DataFrame:
    """Mean and standard deviation across seeds, per endpoint."""
    return (
        per_seed.groupby(["endpoint", "metric"], as_index=False)["value"]
        .agg(["mean", "std", "min", "max"])
        .rename(columns={"mean": "value_mean", "std": "value_std"})
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints", nargs="*", default=list(PK_ENDPOINTS.values()))
    parser.add_argument("--seeds", nargs="*", type=int, default=SEEDS)
    parser.add_argument("--tag", default="baseline_lgbm")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    name_of = {key: name for name, key in PK_ENDPOINTS.items()}

    frames = [run_endpoint(name_of[k], k, args.seeds) for k in args.endpoints]
    per_seed = pd.concat(frames, ignore_index=True)
    summary = summarise(per_seed)

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / f"{args.tag}_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / f"{args.tag}_summary.csv", index=False)
    (RESULTS_ROOT / f"{args.tag}_config.json").write_text(
        json.dumps({"lgbm": LGBM_PARAMS, "seeds": args.seeds}, indent=2)
    )

    print("\n=== Per-seed ===")
    print(per_seed.to_string(index=False))
    print("\n=== Summary (mean +/- sd over seeds) ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
