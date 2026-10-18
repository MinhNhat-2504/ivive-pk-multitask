"""Train on the benchmark, test on drugs approved after it was compiled.

This is the evaluation the whole external-set effort exists for. Models are
fitted on the TDC training data and scored on compounds first approved from
2017 onward, hand-verified against their own label text and disjoint from both
the benchmark and PKSmart's released external sets.

Two endpoints only. `notes/05` established that regulatory documents carry no
intrinsic clearance in assay units at all, and `notes/11` that only 21
compounds yield a true systemic clearance, so the in vitro heads and the
clearance/volume pair cannot be validated this way. Protein binding and
half-life can.

Censoring is handled rather than ignored. 28 of the 129 external binding values
are inequalities -- "greater than 99%" and the like -- and scoring a prediction
against the bound as though it were a measurement would be the same mistake the
benchmark makes. A right-censored point counts as correct when the prediction
lies above the bound, and the uncensored subset is reported alongside.
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

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.data.pk_extract import EXTERNAL_ROOT
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)

#: external file -> (TDC benchmark, log-transform the target?, unit)
ENDPOINTS = {
    "ppb": ("PPBR_AZ", False, "% bound"),
    "half_life": ("Half_Life_Obach", True, "hours"),
}
N_BOOT = 10_000


def load_external(key: str) -> pd.DataFrame:
    frame = pd.read_csv(EXTERNAL_ROOT / f"{key}_final.csv")
    frame["smiles"] = frame["parent_smiles"].map(canonical_smiles)
    frame = frame.dropna(subset=["smiles", "value"]).drop_duplicates("smiles")
    frame["censoring"] = frame["censoring"].fillna("none")
    return frame.reset_index(drop=True)


def censored_absolute_error(
    truth: np.ndarray, prediction: np.ndarray, censoring: np.ndarray
) -> np.ndarray:
    """Absolute error, with no penalty for falling on the right side of a bound.

    A label of "> 99%" is satisfied by any prediction above 99, so the error
    there is zero rather than |prediction - 99|. Treating the bound as a point
    estimate would systematically penalise exactly the compounds the assay
    could not resolve.
    """
    error = np.abs(truth - prediction)
    error = np.where((censoring == "gt") & (prediction >= truth), 0.0, error)
    error = np.where((censoring == "lt") & (prediction <= truth), 0.0, error)
    return error


def evaluate(key: str, seeds: list[int]) -> pd.DataFrame:
    benchmark_name, use_log, unit = ENDPOINTS[key]
    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get(benchmark_name)

    external = load_external(key)
    train_smiles = {canonical_smiles(s) for s in benchmark["train_val"]["Drug"]}
    leaked = external["smiles"].isin(train_smiles)
    if leaked.any():
        log.warning("%s: dropping %d external molecules present in training", key, int(leaked.sum()))
    external = external.loc[~leaked].reset_index(drop=True)

    x_train, ok_train = featurise(benchmark["train_val"]["Drug"].tolist())
    y_train = benchmark["train_val"]["Y"].to_numpy()[ok_train]
    x_bench, ok_bench = featurise(benchmark["test"]["Drug"].tolist())
    y_bench = benchmark["test"]["Y"].to_numpy()[ok_bench]
    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)
    y_ext = external["value"].to_numpy()
    censoring = external["censoring"].to_numpy()
    uncensored = censoring == "none"

    forward = (lambda v: np.log10(np.clip(v, 1e-6, None))) if use_log else (lambda v: v)
    inverse = (lambda v: 10.0**v) if use_log else (lambda v: v)

    rows = []
    for seed in seeds:
        model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
        model.fit(x_train, forward(y_train))
        on_bench = inverse(model.predict(x_bench))
        on_ext = inverse(model.predict(x_ext))

        rows.append(
            {
                "endpoint": key,
                "seed": seed,
                "n_train": int(ok_train.sum()),
                "n_benchmark_test": int(y_bench.size),
                "n_external": int(y_ext.size),
                "n_external_uncensored": int(uncensored.sum()),
                "benchmark_spearman": float(spearmanr(y_bench, on_bench).statistic),
                "benchmark_mae": float(mean_absolute_error(y_bench, on_bench)),
                "external_spearman": float(spearmanr(y_ext, on_ext).statistic),
                "external_mae_censor_aware": float(
                    censored_absolute_error(y_ext, on_ext, censoring).mean()
                ),
                "external_mae_naive": float(mean_absolute_error(y_ext, on_ext)),
                "external_spearman_uncensored": float(
                    spearmanr(y_ext[uncensored], on_ext[uncensored]).statistic
                ),
                "external_mae_uncensored": float(
                    mean_absolute_error(y_ext[uncensored], on_ext[uncensored])
                ),
                # Half-life spans three orders of magnitude, so an error in raw
                # hours is dominated by the long tail. The median ratio of
                # predicted to observed says how many fold out a typical
                # prediction is, which is how the PK literature reports it.
                "external_fold_error": float(
                    np.median(
                        10.0
                        ** np.abs(
                            np.log10(np.clip(on_ext, 1e-6, None))
                            - np.log10(np.clip(y_ext, 1e-6, None))
                        )
                    )
                ),
                "benchmark_fold_error": float(
                    np.median(
                        10.0
                        ** np.abs(
                            np.log10(np.clip(on_bench, 1e-6, None))
                            - np.log10(np.clip(y_bench, 1e-6, None))
                        )
                    )
                ),
            }
        )
        log.info(
            "%s seed %2d: benchmark rho %.3f -> external rho %.3f",
            key, seed, rows[-1]["benchmark_spearman"], rows[-1]["external_spearman"],
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 11)))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    frames = [evaluate(key, args.seeds) for key in ENDPOINTS]
    per_seed = pd.concat(frames, ignore_index=True)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "external_validation_per_seed.csv", index=False)

    rng = np.random.default_rng(0)
    rows = []
    for (endpoint,), group in per_seed.groupby(["endpoint"]):
        row = {"endpoint": endpoint}
        for column in group.columns:
            if column.startswith(("benchmark_", "external_")):
                values = group[column].to_numpy()
                row[column] = values.mean()
                row[f"{column}_sd"] = values.std(ddof=1)
        drop = group["benchmark_spearman"].to_numpy() - group["external_spearman"].to_numpy()
        boot = np.array([rng.choice(drop, drop.size, replace=True).mean() for _ in range(N_BOOT)])
        row["spearman_drop"] = drop.mean()
        row["drop_ci_lo"] = float(np.percentile(boot, 2.5))
        row["drop_ci_hi"] = float(np.percentile(boot, 97.5))
        row["n_external"] = int(group["n_external"].iloc[0])
        row["n_external_uncensored"] = int(group["n_external_uncensored"].iloc[0])
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS_ROOT / "external_validation_summary.csv", index=False)

    for _, row in summary.iterrows():
        unit = ENDPOINTS[row["endpoint"]][2]
        print(f"\n=== {row['endpoint']}  (external n={int(row['n_external'])}, "
              f"uncensored {int(row['n_external_uncensored'])}) ===")
        print(f"  benchmark test   Spearman {row['benchmark_spearman']:.3f} "
              f"± {row['benchmark_spearman_sd']:.3f}   MAE {row['benchmark_mae']:.2f} {unit}")
        print(f"  external         Spearman {row['external_spearman']:.3f} "
              f"± {row['external_spearman_sd']:.3f}   MAE {row['external_mae_censor_aware']:.2f} "
              f"(censor-aware) / {row['external_mae_naive']:.2f} (naive)")
        print(f"  external, uncens Spearman {row['external_spearman_uncensored']:.3f}   "
              f"MAE {row['external_mae_uncensored']:.2f}")
        print(f"  median fold-error  benchmark {row['benchmark_fold_error']:.2f}x  "
              f"-> external {row['external_fold_error']:.2f}x")
        print(f"  DROP in Spearman {row['spearman_drop']:+.3f} "
              f"[{row['drop_ci_lo']:+.3f}, {row['drop_ci_hi']:+.3f}]")


if __name__ == "__main__":
    main()
