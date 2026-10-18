"""All four hyperparameters at once, to rule out interactions.

`notes/09` and `notes/17` swept each knob around the defaults, one at a time.
That cannot see an interaction: a parameter can look irrelevant alone and still
matter in combination. A full factorial over four knobs is 540 configurations,
so the space is sampled instead -- a randomised design covers interactions far
more efficiently than a grid at the same budget, and the question here is
whether any corner of the space behaves differently, not which corner is best.

The defaults are included explicitly so the reported configuration can be
located within the sample.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr

from src.data.chembl_az import DECLARED_RANGE
from src.data.ood_split import butina_clusters, fingerprints
from src.experiments.human_only import build, scaffold_groups
from src.experiments.ood_generalisation import assign_test
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)

DEFAULTS = {"pseudo_weight": 0.3, "n_rounds": 2, "label_weight": 0.05, "min_pairs": 20}
LEVELS = {
    "pseudo_weight": [0.1, 0.2, 0.3, 0.5, 0.8, 1.0],
    "n_rounds": [2, 3, 4],
    "label_weight": [0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
    "min_pairs": [5, 10, 20, 30, 50],
}


def sample_configs(n: int, seed: int = 0) -> list[dict]:
    """A random sample of the joint space, with the defaults included."""
    rng = np.random.default_rng(seed)
    configs = [dict(DEFAULTS)]
    seen = {tuple(DEFAULTS.values())}
    while len(configs) < n:
        candidate = {key: rng.choice(values).item() for key, values in LEVELS.items()}
        key = tuple(candidate.values())
        if key not in seen:
            seen.add(key)
            configs.append(candidate)
    return configs


def run_seed(hep, mic, groups, cluster_of, seed, fractions, configs):
    test_mask = assign_test(groups, len(hep), test_size=0.2, seed=seed)
    train_all, test = hep.loc[~test_mask], hep.loc[test_mask]
    test_clusters = {cluster_of[s] for s in test["smiles"] if s in cluster_of}
    blocked = {s for s, c in cluster_of.items() if c in test_clusters}
    mic_safe = mic[~mic["smiles"].isin(blocked)]

    x_test, ok_test = featurise(test["smiles"].tolist())
    y_test = test["cl_hep"].to_numpy()[ok_test]
    if y_test.size < 20:
        return []

    lo, hi = DECLARED_RANGE["cl_hep_human"]
    mic_lo, mic_hi = DECLARED_RANGE["cl_mic_human"]
    rng = np.random.default_rng(seed)
    rows = []

    for fraction in fractions:
        keep = rng.random(len(train_all)) < fraction
        hep_train = train_all.loc[keep]
        if len(hep_train) < 40:
            continue

        x_a, ok_a = featurise(hep_train["smiles"].tolist())
        y_a = np.log10(np.clip(hep_train["cl_hep"].to_numpy()[ok_a] * HEPATOCYTE_SCALE, 1e-6, None))
        baseline = float(
            spearmanr(
                y_test,
                LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_a, y_a).predict(x_test),
            ).statistic
        )

        pool = pd.merge(hep_train, mic_safe, on="smiles", how="outer")
        x_b, ok_b = featurise(pool["smiles"].tolist())
        pool = pool.loc[ok_b]
        hep_values = pool["cl_hep"].to_numpy(dtype=float)
        mic_values = pool["cl_mic"].to_numpy(dtype=float)
        censored_hep = np.isin(hep_values, [lo, hi])
        censored_mic = np.isin(mic_values, [mic_lo, mic_hi])

        for index, config in enumerate(configs):
            model = ConstrainedMultitask(
                pseudo_weight=config["pseudo_weight"],
                n_rounds=int(config["n_rounds"]),
                label_weight=config["label_weight"],
                min_coupling_pairs=int(config["min_pairs"]),
                random_state=seed,
                **LGBM_PARAMS,
            ).fit(x_b, hep_values, mic_values, censored_hep, censored_mic)
            rows.append(
                {
                    "seed": seed,
                    "fraction": fraction,
                    "config": index,
                    **config,
                    "single_task": baseline,
                    "constrained": float(
                        spearmanr(y_test, model.predict(x_test, "hep")).statistic
                    ),
                }
            )
        log.info("seed %2d frac %.2f done (%d configs)", seed, fraction, len(configs))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 26)))
    parser.add_argument("--fractions", nargs="*", type=float, default=[0.25, 1.0])
    parser.add_argument("--n-configs", type=int, default=24)
    parser.add_argument("--tag", default="joint_sweep")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    configs = sample_configs(args.n_configs)
    hep, mic = build()
    bits, ok = fingerprints(hep["smiles"].tolist())
    hep = hep.loc[ok].reset_index(drop=True)
    groups = scaffold_groups(hep["smiles"].tolist())

    union = sorted(set(hep["smiles"]) | set(mic["smiles"]))
    union_bits, union_ok = fingerprints(union)
    union = [s for s, keep in zip(union, union_ok) if keep]
    cluster_of = {
        union[i]: index
        for index, members in enumerate(butina_clusters(union_bits))
        for i in members
    }

    records = []
    for seed in args.seeds:
        records.extend(run_seed(hep, mic, groups, cluster_of, seed, args.fractions, configs))

    per_seed = pd.DataFrame(records)
    per_seed["delta"] = per_seed["constrained"] - per_seed["single_task"]
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / f"{args.tag}_per_seed.csv", index=False)

    summary = (
        per_seed.groupby(["fraction", "config", "pseudo_weight", "n_rounds",
                          "label_weight", "min_pairs"])["delta"]
        .agg(["mean", "std", lambda s: int((s > 0).sum()), "size"])
        .rename(columns={"<lambda_0>": "wins", "size": "n_seeds"})
        .reset_index()
    )
    summary.to_csv(RESULTS_ROOT / f"{args.tag}_summary.csv", index=False)

    for fraction, group in summary.groupby("fraction"):
        worst = group.nsmallest(3, "mean")
        best = group.nlargest(3, "mean")
        print(f"\n=== {fraction:.0%} labels, {len(group)} sampled configurations ===")
        print(f"  delta: min {group['mean'].min():.4f}  "
              f"median {group['mean'].median():.4f}  max {group['mean'].max():.4f}  "
              f"spread {group['mean'].max() - group['mean'].min():.4f}")
        print(f"  configurations with delta <= 0: "
              f"{int((group['mean'] <= 0).sum())} of {len(group)}")
        print("\n  worst three:")
        print(worst[["pseudo_weight", "n_rounds", "label_weight", "min_pairs", "mean", "wins"]]
              .round(4).to_string(index=False))
        print("  best three:")
        print(best[["pseudo_weight", "n_rounds", "label_weight", "min_pairs", "mean", "wins"]]
              .round(4).to_string(index=False))


if __name__ == "__main__":
    main()
