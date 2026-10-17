"""Does the IVIVE constraint buy accuracy when hepatocyte labels are scarce?

This is the claim `notes/06` left standing. The constraint cannot raise the
leaderboard score -- its own error, 2-3 fold, is the size of the assay noise --
but it can carry information from one task to another across molecules where
the labels do not bridge. If that is real, it shows up as the constrained model
needing fewer labels to reach the same accuracy.

Protocol, per seed and per label fraction:

    arm A  single-task LightGBM on the retained fraction of cl_hep labels
    arm B  constrained model, same cl_hep labels plus the cl_mic labels and
           the learned coupling between the two assay systems

Both are scored on the identical cl_hep test fold, and the difference is taken
within the split, as Section 4 requires.

Leakage control, which this comparison lives or dies by: 682 molecules carry
both labels. A molecule in the cl_hep test fold whose cl_mic label sat in
training would let arm B see a near-measurement of its own target. Every such
molecule is removed from arm B's cl_mic training data, and the count removed is
reported so the control is auditable rather than asserted.
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
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)
FRACTIONS = [0.25, 0.5, 0.75, 1.0]
N_BOOT = 10_000


def load_mic_labels() -> pd.DataFrame:
    frame = ADME(name="Clearance_Microsome_AZ", path=str(DATA_ROOT / "raw")).get_data()
    frame = frame.rename(columns={"Drug": "smiles", "Y": "cl_mic"})[["smiles", "cl_mic"]]
    frame["smiles"] = frame["smiles"].map(canonical_smiles)
    return frame.dropna(subset=["smiles"]).groupby("smiles", as_index=False)["cl_mic"].mean()


def run_seed(
    hep_dataset: ADME, mic_labels: pd.DataFrame, seed: int, fractions: list[float]
) -> list[dict]:
    split = hep_dataset.get_split(method="scaffold", seed=seed)
    train = pd.concat([split["train"], split["valid"]], ignore_index=True)
    test = split["test"]

    train = train.assign(smiles=train["Drug"].map(canonical_smiles)).dropna(subset=["smiles"])
    test = test.assign(smiles=test["Drug"].map(canonical_smiles)).dropna(subset=["smiles"])
    test_smiles = set(test["smiles"])

    # Leakage control: a microsome label on a hepatocyte test molecule is
    # almost a measurement of the target, so it cannot be in training.
    mic_safe = mic_labels[~mic_labels["smiles"].isin(test_smiles)]
    n_dropped = len(mic_labels) - len(mic_safe)

    x_test, ok_test = featurise(test["smiles"].tolist())
    y_test = test["Y"].to_numpy()[ok_test]

    rng = np.random.default_rng(seed)
    rows = []

    for fraction in fractions:
        keep = rng.random(len(train)) < fraction
        hep_train = train.loc[keep, ["smiles", "Y"]].rename(columns={"Y": "cl_hep"})
        if len(hep_train) < 50:
            continue

        # Arm A sees only the retained hepatocyte labels.
        x_a, ok_a = featurise(hep_train["smiles"].tolist())
        y_a = np.log10(np.clip(hep_train["cl_hep"].to_numpy()[ok_a] * HEPATOCYTE_SCALE, 1e-6, None))
        arm_a = LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_a, y_a)
        score_a = float(spearmanr(y_test, arm_a.predict(x_test)).statistic)

        # Arm B additionally sees microsome labels, on the union of molecules.
        pool = pd.merge(hep_train, mic_safe, on="smiles", how="outer")
        x_b, ok_b = featurise(pool["smiles"].tolist())
        pool = pool.loc[ok_b]
        lo, hi = ASSAY_LIMITS["cl_hep"]
        arm_b = ConstrainedMultitask(random_state=seed, **LGBM_PARAMS).fit(
            x_b,
            pool["cl_hep"].to_numpy(dtype=float),
            pool["cl_mic"].to_numpy(dtype=float),
            censored_hep=np.isin(pool["cl_hep"].to_numpy(), [lo, hi]),
            censored_mic=np.isin(pool["cl_mic"].to_numpy(), ASSAY_LIMITS["cl_mic"]),
        )
        score_b = float(spearmanr(y_test, arm_b.predict(x_test, "hep")).statistic)

        rows.append(
            {
                "seed": seed,
                "fraction": fraction,
                "n_hep_labels": int(ok_a.sum()),
                "n_mic_labels": int(np.isfinite(pool["cl_mic"].to_numpy(dtype=float)).sum()),
                "n_mic_dropped_for_leakage": n_dropped,
                "single_task": score_a,
                "constrained": score_b,
            }
        )
        log.info(
            "seed %2d frac %.2f: single %.4f  constrained %.4f  (delta %+.4f)",
            seed, fraction, score_a, score_b, score_b - score_a,
        )
    return rows


def summarise(per_seed: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for fraction, group in per_seed.groupby("fraction"):
        diff = (group["constrained"] - group["single_task"]).to_numpy()
        boot = np.array([rng.choice(diff, diff.size, replace=True).mean() for _ in range(N_BOOT)])
        try:
            p_value = float(wilcoxon(diff).pvalue)
        except ValueError:
            p_value = float("nan")
        rows.append(
            {
                "fraction": fraction,
                "n_hep_labels": int(group["n_hep_labels"].mean()),
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
    per_seed.to_csv(RESULTS_ROOT / "label_efficiency_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / "label_efficiency_summary.csv", index=False)

    print("\n=== Constrained vs single-task on cl_hep, paired per split ===")
    print(summary.round(4).to_string(index=False))
    print(
        f"\nMicrosome labels withheld for leakage control: "
        f"{int(per_seed['n_mic_dropped_for_leakage'].mean())} per seed on average."
    )


if __name__ == "__main__":
    main()
