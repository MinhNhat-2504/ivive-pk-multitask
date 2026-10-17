"""The label-efficiency experiment again, on human data only.

Section 7 trained on TDC's `Clearance_Hepatocyte_AZ`, which turned out to be
about two thirds rat. The paired comparisons survive that -- every arm saw the
same data, so the contrast is unaffected -- but the pharmacological reading
does not: the learned coupling was bridging human microsomes to a human/rat
mixture, and part of what it absorbed was a species difference rather than a
difference between in vitro systems.

So the experiment is repeated on the deposition's human hepatocyte assay alone
(408 compounds, CHEMBL3301372) against human microsomes (1102,
CHEMBL3301370). It is a much smaller and harder setting -- 187 molecules carry
both labels against the 682 the pooled set appeared to offer -- and that is the
point: if the constraint still pays here, it pays for the reason claimed.

Splits are Bemis-Murcko scaffold groups, matching the setting of Section 7 so
the two are comparable. A cluster split was tried first and is too hard to
measure anything on at this size: single-task Spearman falls to 0.047 +/- 0.161
against 0.374 +/- 0.081 under a random split, so the model has no measurable
ability to extrapolate to unseen chemical series from 408 compounds and no
constraint effect could be resolved against that noise. That is reported as a
result in its own right rather than worked around.

Microsome labels are still withheld at cluster level, the strict control from
`notes/10`, since that is about leakage rather than split difficulty.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr, wilcoxon

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from src.data.chembl_az import DECLARED_RANGE, load
from src.data.ood_split import butina_clusters, fingerprints
from src.experiments.ood_generalisation import assign_test
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.constrained import ConstrainedMultitask
from src.models.ivive import HEPATOCYTE_SCALE

log = logging.getLogger(__name__)
N_BOOT = 10_000
FRACTIONS = [0.25, 0.5, 0.75, 1.0]


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Human hepatocyte targets and human microsome sources, deduplicated."""
    frame = load()
    hep = (
        frame[frame["assay"] == "cl_hep_human"]
        .groupby("smiles", as_index=False)["value"]
        .mean()
        .rename(columns={"value": "cl_hep"})
    )
    mic = (
        frame[frame["assay"] == "cl_mic_human"]
        .groupby("smiles", as_index=False)["value"]
        .mean()
        .rename(columns={"value": "cl_mic"})
    )
    return hep.reset_index(drop=True), mic


def scaffold_groups(smiles: list[str]) -> list[list[int]]:
    """Group molecule indices by Bemis-Murcko framework."""
    groups: dict[str, list[int]] = {}
    for index, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        try:
            key = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) if mol else ""
        except Exception:
            key = ""
        groups.setdefault(key or f"__singleton_{index}", []).append(index)
    return list(groups.values())


def run_seed(
    hep: pd.DataFrame,
    mic: pd.DataFrame,
    clusters: list[list[int]],
    cluster_of: dict[str, int],
    seed: int,
    fractions: list[float],
) -> list[dict]:
    test_mask = assign_test(clusters, len(hep), test_size=0.2, seed=seed)
    train_all, test = hep.loc[~test_mask], hep.loc[test_mask]

    # Strict control: withhold every microsome molecule sharing a cluster with
    # the test set, not merely the test molecules themselves.
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
        scores = {
            "single_task": float(
                spearmanr(
                    y_test,
                    LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_a, y_a).predict(x_test),
                ).statistic
            )
        }

        pool = pd.merge(hep_train, mic_safe, on="smiles", how="outer")
        x_b, ok_b = featurise(pool["smiles"].tolist())
        pool = pool.loc[ok_b]
        hep_values = pool["cl_hep"].to_numpy(dtype=float)
        mic_values = pool["cl_mic"].to_numpy(dtype=float)

        for arm, mode in (("constrained", "deming"), ("constrained_identity", "identity")):
            model = ConstrainedMultitask(
                coupling_mode=mode, random_state=seed, **LGBM_PARAMS
            ).fit(
                x_b,
                hep_values,
                mic_values,
                censored_hep=np.isin(hep_values, [lo, hi]),
                censored_mic=np.isin(mic_values, [mic_lo, mic_hi]),
            )
            scores[arm] = float(spearmanr(y_test, model.predict(x_test, "hep")).statistic)

        rows.append(
            {
                "seed": seed,
                "fraction": fraction,
                "n_hep_labels": int(ok_a.sum()),
                "n_mic_labels": int(np.isfinite(mic_values).sum()),
                "n_test": int(y_test.size),
                **scores,
            }
        )
        log.info(
            "seed %2d frac %.2f (n=%3d): single %.4f  deming %.4f  identity %.4f",
            seed, fraction, int(ok_a.sum()),
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
                    "n_hep_labels": int(group["n_hep_labels"].mean()),
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
    parser.add_argument("--split", choices=["scaffold", "cluster"], default="scaffold")
    parser.add_argument("--tag", default="human_only")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    hep, mic = build()
    log.info("human hepatocyte %d, human microsome %d, overlap %d",
             len(hep), len(mic), len(set(hep["smiles"]) & set(mic["smiles"])))

    bits, ok = fingerprints(hep["smiles"].tolist())
    hep = hep.loc[ok].reset_index(drop=True)
    groups = (
        scaffold_groups(hep["smiles"].tolist())
        if args.split == "scaffold"
        else butina_clusters(bits)
    )
    log.info("%s split: %d groups over %d molecules", args.split, len(groups), len(hep))

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

    print(f"\nmicrosome labels available after cluster-level withholding: "
          f"{int(per_seed['n_mic_labels'].mean())} of {len(mic)}")
    for fraction, group in summary.groupby("fraction"):
        print(f"\n=== human only, {fraction:.0%} labels (n={int(group['n_hep_labels'].iloc[0])}) ===")
        print(
            group[["comparison", "arm_mean", "reference_mean", "mean_delta", "ci_lo", "ci_hi", "wins", "wilcoxon_p"]]
            .round(4)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
