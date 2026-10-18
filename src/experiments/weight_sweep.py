"""The last two unswept knobs, on the human-only setting.

`notes/09` swept `pseudo_weight` and the number of rounds. Two remained, and
both turned out to matter more than their obscurity suggested:

* **the projection weights.** The constraint pulls a molecule's two clearance
  predictions onto the constraint line, weighted so a head holding a
  measurement barely moves while a purely predicted one absorbs the correction.
  The ratio was fixed at 0.05 against 1.0 and justified by a within-molecule
  spread that `notes/15` retracted, so it now needs a sweep rather than a
  citation. A label weight of 1.0 removes the distinction entirely.

* **the fallback threshold.** Below a minimum number of molecules carrying both
  labels, the coupling is not fitted and the model reverts to pure
  physiological scaling. `notes/16` showed this rule, fixed at 20, is what the
  25%-label ablation row actually measures: the paired count there is 15.9 on
  average, under the threshold in 82% of runs, so both arms become the same
  model. Sweeping it separates the rule from the result.

Run on human-only data, which `notes/16` makes the primary setting, paired per
split against the same single-task baseline.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr, wilcoxon

from src.data.chembl_az import DECLARED_RANGE
from src.data.ood_split import butina_clusters, fingerprints
from src.experiments.human_only import build, scaffold_groups
from src.experiments.ood_generalisation import assign_test
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)
N_BOOT = 10_000

#: (label_weight, min_coupling_pairs). The defaults sit mid-range in both.
GRID = [
    (0.01, 20), (0.05, 20), (0.10, 20), (0.25, 20), (0.50, 20), (1.00, 20),
    (0.05, 5), (0.05, 10), (0.05, 30), (0.05, 50),
]
FRACTIONS = [0.25, 1.0]


def run_seed(hep, mic, groups, cluster_of, seed, fractions):
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

        for label_weight, min_pairs in GRID:
            model = ConstrainedMultitask(
                label_weight=label_weight,
                min_coupling_pairs=min_pairs,
                random_state=seed,
                **LGBM_PARAMS,
            ).fit(x_b, hep_values, mic_values, censored_hep, censored_mic)
            rows.append(
                {
                    "seed": seed,
                    "fraction": fraction,
                    "label_weight": label_weight,
                    "min_pairs": min_pairs,
                    "n_hep_labels": int(ok_a.sum()),
                    "coupling_n": model.coupling_.n_fitted,
                    "coupling_slope": model.coupling_.slope,
                    "single_task": baseline,
                    "constrained": float(
                        spearmanr(y_test, model.predict(x_test, "hep")).statistic
                    ),
                }
            )
        log.info("seed %2d frac %.2f done (%d configs)", seed, fraction, len(GRID))
    return rows


def summarise(per_seed: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for (fraction, weight, pairs), group in per_seed.groupby(
        ["fraction", "label_weight", "min_pairs"]
    ):
        diff = (group["constrained"] - group["single_task"]).to_numpy()
        boot = np.array([rng.choice(diff, diff.size, replace=True).mean() for _ in range(N_BOOT)])
        try:
            p_value = float(wilcoxon(diff).pvalue)
        except ValueError:
            p_value = float("nan")
        rows.append(
            {
                "fraction": fraction,
                "label_weight": weight,
                "min_pairs": pairs,
                "constrained": group["constrained"].mean(),
                "mean_delta": diff.mean(),
                "ci_lo": float(np.percentile(boot, 2.5)),
                "ci_hi": float(np.percentile(boot, 97.5)),
                "wins": int((diff > 0).sum()),
                "n_seeds": diff.size,
                "coupling_fitted_pct": 100.0 * (group["coupling_slope"] != 1.0).mean(),
                "wilcoxon_p": p_value,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 26)))
    parser.add_argument("--fractions", nargs="*", type=float, default=FRACTIONS)
    parser.add_argument("--tag", default="weight_sweep")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

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
        records.extend(run_seed(hep, mic, groups, cluster_of, seed, args.fractions))

    per_seed = pd.DataFrame(records)
    summary = summarise(per_seed)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / f"{args.tag}_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / f"{args.tag}_summary.csv", index=False)

    for fraction, group in summary.groupby("fraction"):
        print(f"\n=== human only, {fraction:.0%} labels — constrained minus single-task ===")
        print(
            group[["label_weight", "min_pairs", "mean_delta", "ci_lo", "ci_hi",
                   "wins", "coupling_fitted_pct", "wilcoxon_p"]]
            .round(4)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
