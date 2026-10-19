"""Does the target scale still matter once both scales are tuned?

Every comparison so far used one untuned configuration, which invites the
obvious objection: the percent-bound target may simply be under-served by
hyperparameters chosen without it in mind, and a tuned percent model might close
the gap.

The test has to be fair in two respects. Both scales get the same search budget
over the same grid. And both are selected on the benchmark's own metric --
validation mean absolute error on percent bound -- because that is what someone
submitting to the leaderboard would tune for. Selecting on the unbound-fraction
error instead would hand the transformed targets the comparison by construction.

Predictions, recorded before running:

    1.  Tuning improves both scales on mean absolute error, by a similar small
        amount, because the search is optimising exactly that quantity.
    2.  Tuning does not close the fu gap. The selection criterion is blind to
        the region where fu varies -- that is the paper's argument -- so no
        amount of searching under it can recover what the target scale throws
        away.
    3.  A tuned percent model may therefore end up with a *better* mean absolute
        error than an untuned logit model while still being far worse as a
        source of fu. If that happens it is the cleanest statement of the
        problem the paper describes.

If prediction 2 fails -- if a tuned percent model reaches the fu accuracy of a
logit model -- then the scale was a proxy for hyperparameters and the paper's
central claim is wrong.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.experiments.target_scale import HIGH_BINDING, SCALES, to_fu
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)

SEEDS = [1, 2, 3, 4, 5]
N_CONFIGS = 16

#: The search space. Deliberately centred on the untuned configuration rather
#: than around it, so that "tuned" cannot mean "moved somewhere the baseline
#: was never allowed to go".
SEARCH_SPACE = {
    "n_estimators": [300, 600, 1200],
    "learning_rate": [0.01, 0.03, 0.06],
    "num_leaves": [15, 31, 63],
    "min_child_samples": [5, 10, 25],
    "colsample_bytree": [0.4, 0.6, 0.8],
    "subsample": [0.6, 0.8, 1.0],
    "reg_lambda": [0.1, 1.0, 10.0],
}


def sample_configs(rng: np.random.Generator) -> list[dict]:
    """Random configurations, with the untuned one always included."""
    configs = [dict(LGBM_PARAMS)]
    while len(configs) < N_CONFIGS:
        candidate = dict(LGBM_PARAMS)
        for key, values in SEARCH_SPACE.items():
            candidate[key] = values[int(rng.integers(len(values)))]
        candidate["subsample_freq"] = 1
        configs.append(candidate)
    return configs


def score(truth: np.ndarray, pred: np.ndarray) -> dict:
    fold = np.maximum(to_fu(pred) / to_fu(truth), to_fu(truth) / to_fu(pred))
    high = truth >= HIGH_BINDING
    return {
        "mae": float(np.mean(np.abs(pred - truth))),
        "spearman": float(spearmanr(truth, pred).statistic),
        "fu_fold_median_high": float(np.median(fold[high])),
        "frac_within_2fold_high": float(np.mean(fold[high] <= 2.0)),
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    external = load_external("ppb")
    train_smiles = {canonical_smiles(s) for s in benchmark["train_val"]["Drug"]}
    external = external[~external["smiles"].isin(train_smiles)].reset_index(drop=True)

    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)
    y_ext = external["value"].to_numpy()
    x_test, ok_test = featurise(benchmark["test"]["Drug"].tolist())
    y_test = benchmark["test"]["Y"].to_numpy()[ok_test]

    rows = []
    for seed in SEEDS:
        train, valid = group.get_train_valid_split(
            benchmark="PPBR_AZ", split_type="default", seed=seed
        )
        x_train, ok_train = featurise(train["Drug"].tolist())
        y_train = train["Y"].to_numpy()[ok_train]
        x_valid, ok_valid = featurise(valid["Drug"].tolist())
        y_valid = valid["Y"].to_numpy()[ok_valid]

        configs = sample_configs(np.random.default_rng(seed))
        for scale, (forward, inverse) in SCALES.items():
            best, best_mae = None, np.inf
            for index, config in enumerate(configs):
                model = LGBMRegressor(random_state=seed, **config)
                model.fit(x_train, forward(y_train))
                # Selected on the leaderboard's metric, on held-out validation.
                mae = float(np.mean(np.abs(inverse(model.predict(x_valid)) - y_valid)))
                if mae < best_mae:
                    best, best_mae, best_index = model, mae, index

            rows.append({
                "seed": seed, "scale": scale, "tuning": "tuned",
                "config_index": best_index, "valid_mae": best_mae,
                **{f"test_{k}": v for k, v in
                   score(y_test, inverse(best.predict(x_test))).items()},
                **{f"ext_{k}": v for k, v in
                   score(y_ext, inverse(best.predict(x_ext))).items()},
            })

            baseline = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
            baseline.fit(x_train, forward(y_train))
            rows.append({
                "seed": seed, "scale": scale, "tuning": "untuned",
                "config_index": 0,
                "valid_mae": float(np.mean(
                    np.abs(inverse(baseline.predict(x_valid)) - y_valid))),
                **{f"test_{k}": v for k, v in
                   score(y_test, inverse(baseline.predict(x_test))).items()},
                **{f"ext_{k}": v for k, v in
                   score(y_ext, inverse(baseline.predict(x_ext))).items()},
            })
            log.info("seed %d %-8s tuned config %d, valid MAE %.3f",
                     seed, scale, best_index, best_mae)

    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(RESULTS_ROOT / "tuned_scale_per_seed.csv", index=False)

    print(f"\n=== {N_CONFIGS} configurations, selected on validation MAE, "
          f"{len(SEEDS)} seeds ===")
    for metric in ("test_mae", "test_spearman", "test_fu_fold_median_high",
                   "ext_fu_fold_median_high", "ext_frac_within_2fold_high"):
        table = per_seed.pivot_table(index="scale", columns="tuning",
                                     values=metric, aggfunc="mean")
        print(f"\n{metric}")
        print(table.reindex(list(SCALES))[["untuned", "tuned"]].round(4).to_string())

    print("\n=== The question the paper turns on ===")
    tuned = per_seed[per_seed["tuning"] == "tuned"].groupby("scale")
    untuned = per_seed[per_seed["tuning"] == "untuned"].groupby("scale")
    print(f"  tuned percent    test MAE {tuned.get_group('percent').test_mae.mean():.3f}"
          f"   fu fold {tuned.get_group('percent').ext_fu_fold_median_high.mean():.2f}")
    print(f"  untuned logit    test MAE {untuned.get_group('logit').test_mae.mean():.3f}"
          f"   fu fold {untuned.get_group('logit').ext_fu_fold_median_high.mean():.2f}")
    print(f"  tuned logit      test MAE {tuned.get_group('logit').test_mae.mean():.3f}"
          f"   fu fold {tuned.get_group('logit').ext_fu_fold_median_high.mean():.2f}")


if __name__ == "__main__":
    main()
