"""Does the constraint help where the model has to extrapolate?

Section 7 found the constraint pays when labels are scarce and costs a little
when they are plentiful. That was measured on the benchmark's own scaffold
split, which turns out to be a weak extrapolation test: a third to a half of
its test molecules have a training neighbour above 0.5 Tanimoto.

Here the same three arms are run on cluster splits, where whole chemical series
are held out and almost no test molecule has a close training relative. The
claim under test is the one Section 3 left open -- that a mechanistic
constraint keeps predictions admissible where the data runs out -- and the
sharpest case is full labels, where in-domain the constraint was slightly
harmful.
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
from src.data.ood_split import butina_clusters, describe_split, fingerprints
from src.eval.censoring import ASSAY_LIMITS
from src.experiments.label_efficiency import load_mic_labels
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)
N_BOOT = 10_000
FRACTIONS = [0.5, 1.0]


def assign_test(clusters: list[list[int]], n_total: int, test_size: float, seed: int) -> np.ndarray:
    """Hold out whole clusters until the test set reaches the target size."""
    rng = np.random.default_rng(seed)
    picked: list[int] = []
    target = int(round(test_size * n_total))
    for c in rng.permutation(len(clusters)):
        if len(picked) >= target:
            break
        picked.extend(clusters[c])
    mask = np.zeros(n_total, dtype=bool)
    mask[np.array(picked, dtype=int)] = True
    return mask


def run_seed(
    frame: pd.DataFrame,
    clusters: list[list[int]],
    mic_labels: pd.DataFrame,
    seed: int,
    fractions: list[float],
    mic_cluster_of: dict[str, int] | None = None,
) -> list[dict]:
    test_mask = assign_test(clusters, len(frame), test_size=0.2, seed=seed)
    train_all = frame.loc[~test_mask]
    test = frame.loc[test_mask]

    difficulty = describe_split(
        frame["smiles"].tolist(), ~test_mask, test_mask
    )

    # Excluding only the test molecules themselves is too weak for an
    # extrapolation claim: a microsome label on a close analogue of a test
    # molecule tells the constrained arm about that chemical series, which is
    # the very thing being held out. Under `mic_cluster_of` the exclusion
    # extends to every microsome molecule sharing a cluster with the test set.
    excluded = set(test["smiles"])
    if mic_cluster_of is not None:
        test_clusters = {mic_cluster_of[s] for s in test["smiles"] if s in mic_cluster_of}
        excluded |= {s for s, c in mic_cluster_of.items() if c in test_clusters}
    mic_safe = mic_labels[~mic_labels["smiles"].isin(excluded)]
    x_test, ok_test = featurise(test["smiles"].tolist())
    y_test = test["cl_hep"].to_numpy()[ok_test]

    rng = np.random.default_rng(seed)
    lo, hi = ASSAY_LIMITS["cl_hep"]
    rows = []

    for fraction in fractions:
        keep = rng.random(len(train_all)) < fraction
        hep_train = train_all.loc[keep]
        if len(hep_train) < 50:
            continue

        x_a, ok_a = featurise(hep_train["smiles"].tolist())
        y_a = np.log10(np.clip(hep_train["cl_hep"].to_numpy()[ok_a] * HEPATOCYTE_SCALE, 1e-6, None))
        scores = {
            "single_task": float(
                spearmanr(
                    y_test,
                    LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_a, y_a).predict(x_test),
                ).statistic
            )
        }

        pool = pd.merge(hep_train[["smiles", "cl_hep"]], mic_safe, on="smiles", how="outer")
        x_b, ok_b = featurise(pool["smiles"].tolist())
        pool = pool.loc[ok_b]
        hep_vals = pool["cl_hep"].to_numpy(dtype=float)
        mic_vals = pool["cl_mic"].to_numpy(dtype=float)

        for arm, mode in (("constrained", "deming"), ("constrained_identity", "identity")):
            model = ConstrainedMultitask(
                coupling_mode=mode, random_state=seed, **LGBM_PARAMS
            ).fit(
                x_b,
                hep_vals,
                mic_vals,
                censored_hep=np.isin(hep_vals, [lo, hi]),
                censored_mic=np.isin(mic_vals, ASSAY_LIMITS["cl_mic"]),
            )
            scores[arm] = float(spearmanr(y_test, model.predict(x_test, "hep")).statistic)

        rows.append(
            {
                "seed": seed,
                "fraction": fraction,
                "n_train": int(ok_a.sum()),
                "n_test": int(ok_test.sum()),
                "n_mic_available": int(len(mic_safe)),
                "nn_sim_median": difficulty.get("nn_sim_median", np.nan),
                **scores,
            }
        )
        log.info(
            "seed %2d frac %.2f (nn_sim %.3f): single %.4f  deming %.4f  identity %.4f",
            seed, fraction, difficulty.get("nn_sim_median", np.nan),
            scores["single_task"], scores["constrained"], scores["constrained_identity"],
        )
    return rows


def summarise(per_seed: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for fraction, group in per_seed.groupby("fraction"):
        for arm, reference in (
            ("constrained", "single_task"),
            ("constrained_identity", "single_task"),
            ("constrained", "constrained_identity"),
        ):
            diff = (group[arm] - group[reference]).to_numpy()
            boot = np.array(
                [rng.choice(diff, diff.size, replace=True).mean() for _ in range(N_BOOT)]
            )
            try:
                p_value = float(wilcoxon(diff).pvalue)
            except ValueError:
                p_value = float("nan")
            rows.append(
                {
                    "fraction": fraction,
                    "comparison": f"{arm} - {reference}",
                    "arm_mean": group[arm].mean(),
                    "reference_mean": group[reference].mean(),
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
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also withhold microsome labels sharing a cluster with the test set",
    )
    parser.add_argument("--tag", default="ood")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    raw = ADME(name="Clearance_Hepatocyte_AZ", path=str(DATA_ROOT / "raw")).get_data()
    frame = (
        raw.assign(smiles=raw["Drug"].map(canonical_smiles))
        .dropna(subset=["smiles"])
        .groupby("smiles", as_index=False)["Y"]
        .mean()
        .rename(columns={"Y": "cl_hep"})
        .reset_index(drop=True)
    )

    # Clustering is deterministic, so it is computed once and only the
    # assignment of clusters to the test side varies with the seed.
    log.info("clustering %d molecules", len(frame))
    bits, ok = fingerprints(frame["smiles"].tolist())
    frame = frame.loc[ok].reset_index(drop=True)
    clusters = butina_clusters(bits)
    log.info("%d clusters; largest %d", len(clusters), max(len(c) for c in clusters))

    mic_labels = load_mic_labels()

    mic_cluster_of = None
    if args.strict:
        # Cluster the union so microsome-only molecules carry a cluster label
        # comparable with the hepatocyte molecules the test set is drawn from.
        union = sorted(set(frame["smiles"]) | set(mic_labels["smiles"]))
        log.info("strict mode: clustering union of %d molecules", len(union))
        union_bits, union_ok = fingerprints(union)
        union = [s for s, keep in zip(union, union_ok) if keep]
        mic_cluster_of = {
            union[i]: cluster_index
            for cluster_index, members in enumerate(butina_clusters(union_bits))
            for i in members
        }

    records = []
    for seed in args.seeds:
        records.extend(
            run_seed(frame, clusters, mic_labels, seed, args.fractions, mic_cluster_of)
        )

    per_seed = pd.DataFrame(records)
    summary = summarise(per_seed)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / f"{args.tag}_per_seed.csv", index=False)
    summary.to_csv(RESULTS_ROOT / f"{args.tag}_summary.csv", index=False)

    print(f"\nmedian nearest-neighbour similarity, test to train: "
          f"{per_seed['nn_sim_median'].mean():.3f}")
    for fraction, group in summary.groupby("fraction"):
        print(f"\n=== cluster split, {fraction:.0%} labels ===")
        print(
            group[["comparison", "arm_mean", "reference_mean", "mean_delta", "ci_lo", "ci_hi", "wins", "wilcoxon_p"]]
            .round(4)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
