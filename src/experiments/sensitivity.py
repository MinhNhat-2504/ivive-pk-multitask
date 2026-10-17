"""Is the label-efficiency result an artefact of two untuned knobs?

`pseudo_weight` and `n_rounds` were set once, to 0.3 and 2, and never revisited.
If the gain in Section 7 only exists near those values it is a tuning artefact
and the section has to be rewritten; if it holds across the range, the choice
was incidental. Either answer is worth more than the result taken on trust.

Both knobs are swept against the same single-task baseline, paired per split,
at the two label fractions where Section 7 found an effect.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr, wilcoxon
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.eval.censoring import ASSAY_LIMITS
from src.experiments.label_efficiency import load_mic_labels
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)

#: (pseudo_weight, n_rounds). The default sits in the middle of the sweep so a
#: peak exactly at the default would be visible as such.
GRID = [
    (0.1, 2), (0.2, 2), (0.3, 2), (0.5, 2), (0.8, 2), (1.0, 2),
    (0.3, 3), (0.3, 4),
]
FRACTIONS = [0.25, 0.5]
N_BOOT = 10_000


def run_seed(hep_dataset: ADME, mic_labels: pd.DataFrame, seed: int, fractions: list[float]):
    split = hep_dataset.get_split(method="scaffold", seed=seed)
    train = pd.concat([split["train"], split["valid"]], ignore_index=True)
    test = split["test"]
    train = train.assign(smiles=train["Drug"].map(canonical_smiles)).dropna(subset=["smiles"])
    test = test.assign(smiles=test["Drug"].map(canonical_smiles)).dropna(subset=["smiles"])

    mic_safe = mic_labels[~mic_labels["smiles"].isin(set(test["smiles"]))]
    x_test, ok_test = featurise(test["smiles"].tolist())
    y_test = test["Y"].to_numpy()[ok_test]

    rng = np.random.default_rng(seed)
    rows = []
    lo, hi = ASSAY_LIMITS["cl_hep"]

    for fraction in fractions:
        keep = rng.random(len(train)) < fraction
        hep_train = train.loc[keep, ["smiles", "Y"]].rename(columns={"Y": "cl_hep"})
        if len(hep_train) < 50:
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
        hep_vals = pool["cl_hep"].to_numpy(dtype=float)
        mic_vals = pool["cl_mic"].to_numpy(dtype=float)
        censored_hep = np.isin(hep_vals, [lo, hi])
        censored_mic = np.isin(mic_vals, ASSAY_LIMITS["cl_mic"])

        for pseudo_weight, n_rounds in GRID:
            model = ConstrainedMultitask(
                n_rounds=n_rounds,
                pseudo_weight=pseudo_weight,
                random_state=seed,
                **LGBM_PARAMS,
            ).fit(x_b, hep_vals, mic_vals, censored_hep, censored_mic)
            rows.append(
                {
                    "seed": seed,
                    "fraction": fraction,
                    "pseudo_weight": pseudo_weight,
                    "n_rounds": n_rounds,
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
    for (fraction, weight, rounds), group in per_seed.groupby(
        ["fraction", "pseudo_weight", "n_rounds"]
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
                "pseudo_weight": weight,
                "n_rounds": rounds,
                "single_task": group["single_task"].mean(),
                "constrained": group["constrained"].mean(),
                "mean_delta": diff.mean(),
                "ci_lo": float(np.percentile(boot, 2.5)),
                "ci_hi": float(np.percentile(boot, 97.5)),
                "wins": int((diff > 0).sum()),
                "n_seeds": diff.size,
                "wilcoxon_p": p_value,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", nargs="*", type=int, default=list(range(1, 26)))
    parser.add_argument("--fractions", nargs="*", type=float, default=FRACTIONS)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    hep_dataset = ADME(name="Clearance_Hepatocyte_AZ", path=str(DATA_ROOT / "raw"))
    mic_labels = load_mic_labels()

    records = []
    for seed in args.seeds:
        records.extend(run_seed(hep_dataset, mic_labels, seed, args.fractions))

    per_seed = pd.DataFrame(records)
    summary = summarise(per_seed)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "sensitivity_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / "sensitivity_summary.csv", index=False)

    for fraction, group in summary.groupby("fraction"):
        print(f"\n=== label fraction {fraction:.0%} — constrained minus single-task ===")
        print(
            group[
                ["pseudo_weight", "n_rounds", "constrained", "mean_delta", "ci_lo", "ci_hi", "wins", "wilcoxon_p"]
            ]
            .round(4)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
