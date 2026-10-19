"""Does the target-scale effect survive when the training data is one species?

`PPBR_AZ` turns out to pool five plasma protein binding assays from its source
deposition -- human, rat, dog, mouse and guinea pig -- and 42.4% of its rows
carry a non-human value. Every earlier result in this line of work trained on
that pooled label and was evaluated against human values from approved labelling,
so two explanations were confounded:

    scale      the loss is stated on percent bound while the consumer divides
               by fu, so the optimiser is directed away from the top percentage
               point
    species    the model is partly fitted to rat, dog, mouse and guinea pig
               plasma and asked to predict human plasma

Both would compress predictions toward the interior and both would inflate the
fu error for highly bound compounds. They are separated here by refitting on the
human assay alone, CHEMBL3301365, 1614 records with one record per molecule and
no duplication.

The predictions, recorded before running:

    1.  The scale effect survives. It is an argument about the loss, and the
        loss does not know what species produced the label, so restricting to
        human should not remove it.
    2.  The absolute fu error falls on every scale, because the species mismatch
        was real error and removing it removes error.
    3.  Therefore the previously reported 8.3-fold figure for `percent` is partly
        species mismatch, and the honest attribution needs both numbers.

If prediction 1 fails -- if the scale advantage disappears on human-only data --
then the compression was species pooling all along and the scale argument was an
artefact of a defect in the benchmark.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr, wilcoxon
from tdc.benchmark_group import admet_group

from src.data.chembl_az import load
from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.experiments.target_scale import HIGH_BINDING, SCALES, to_fu
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)

SEEDS = list(range(1, 26))
#: Species assays in the deposition, by ChEMBL id; only the first is human.
HUMAN_ASSAY = "CHEMBL3301365"


def score(truth: np.ndarray, pred: np.ndarray) -> dict:
    fold = np.maximum(to_fu(pred) / to_fu(truth), to_fu(truth) / to_fu(pred))
    high = truth >= HIGH_BINDING
    return {
        "mae_percent": float(np.mean(np.abs(pred - truth))),
        "spearman": float(spearmanr(truth, pred).statistic),
        "fu_fold_median": float(np.median(fold)),
        "fu_fold_median_high": float(np.median(fold[high])),
        "frac_within_2fold_high": float(np.mean(fold[high] <= 2.0)),
        "pred_range": float(pred.max() - pred.min()),
        "n_high": int(high.sum()),
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    external = load_external("ppb")
    upstream = pd.read_csv(RESULTS_ROOT.parent / "results" / "ppb_all_species_upstream.csv")

    # Training sets: the pooled benchmark as distributed, and the human assay alone.
    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    pooled = pd.DataFrame({
        "smiles": [canonical_smiles(s) for s in benchmark["train_val"]["Drug"]],
        "value": benchmark["train_val"]["Y"].to_numpy(),
    })

    az = load()
    human = az[az["assay"] == "ppb_human"][["smiles", "value"]].copy()
    human["smiles"] = human["smiles"].map(canonical_smiles)
    human = human.dropna(subset=["smiles", "value"]).drop_duplicates("smiles")

    # The external set must not appear in either training set.
    external = external[~external["smiles"].isin(set(pooled["smiles"]))]
    external = external[~external["smiles"].isin(set(human["smiles"]))].reset_index(drop=True)
    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)
    y_ext = external["value"].to_numpy()

    # Held-out human rows of the benchmark test partition, for a second view
    # that is not the FDA set.
    test = benchmark["test"].copy()
    test["c"] = [canonical_smiles(s) for s in test["Drug"]]
    human_keys = set(
        zip(upstream.loc[upstream.species == "human", "c"],
            upstream.loc[upstream.species == "human", "value"].round(4))
    )
    nonhuman_keys = set(
        zip(upstream.loc[upstream.species != "human", "c"],
            upstream.loc[upstream.species != "human", "value"].round(4))
    )
    test["k"] = list(zip(test["c"], test["Y"].round(4)))
    test_human = test[test["k"].isin(human_keys) & ~test["k"].isin(nonhuman_keys)]
    test_human = test_human.drop_duplicates("c").reset_index(drop=True)
    x_th, ok_th = featurise(test_human["Drug"].tolist())
    test_human = test_human.loc[ok_th].reset_index(drop=True)
    y_th = test_human["Y"].to_numpy()

    log.info("external %d (%d high); benchmark test, human-only rows %d (%d high)",
             len(y_ext), int((y_ext >= HIGH_BINDING).sum()),
             len(y_th), int((y_th >= HIGH_BINDING).sum()))

    training = {}
    for name, frame in (("pooled (as distributed)", pooled), ("human only", human)):
        # Molecules of the human-only test view are removed from both training
        # sets, so the two are scored on the same held-out molecules.
        frame = frame[~frame["smiles"].isin(set(test_human["c"]))]
        x, ok = featurise(frame["smiles"].tolist())
        training[name] = (x, frame["value"].to_numpy()[ok])
        log.info("%-24s n = %d", name, len(training[name][1]))

    rows = []
    for seed in SEEDS:
        for train_name, (x_train, y_train) in training.items():
            for scale, (forward, inverse) in SCALES.items():
                model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
                model.fit(x_train, forward(y_train))
                for view, (x_eval, y_eval) in (
                    ("external", (x_ext, y_ext)), ("test_human", (x_th, y_th))
                ):
                    rows.append({
                        "seed": seed, "training": train_name,
                        "scale": scale, "view": view,
                        **score(y_eval, inverse(model.predict(x_eval))),
                    })

    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(RESULTS_ROOT / "human_only_scale_per_seed.csv", index=False)

    for view in ("external", "test_human"):
        block = per_seed[per_seed["view"] == view]
        print(f"\n=== {view} ===")
        for metric in ("mae_percent", "spearman", "fu_fold_median_high",
                       "frac_within_2fold_high", "pred_range"):
            table = block.pivot_table(
                index="training", columns="scale", values=metric, aggfunc="mean"
            )[list(SCALES)]
            print(f"\n{metric}")
            print(table.round(3).to_string())

    print("\n=== Scale effect within each training set, paired by seed ===")
    print("(external view, median fu fold error for compounds bound above 99%)")
    block = per_seed[per_seed["view"] == "external"]
    wide = block.pivot_table(index=["training", "seed"], columns="scale",
                             values="fu_fold_median_high")
    for train_name in training:
        sub = wide.loc[train_name]
        for scale in ("log_fu", "logit"):
            delta = (sub[scale] - sub["percent"]).to_numpy()
            print(f"  {train_name:>24s}  {scale:>7s} - percent: {delta.mean():+7.3f}  "
                  f"p={wilcoxon(delta).pvalue:.2g}  better in "
                  f"{int((delta < 0).sum())}/{len(delta)}")

    print("\n=== Species effect within each scale, paired by seed ===")
    wide2 = block.pivot_table(index=["scale", "seed"], columns="training",
                              values="fu_fold_median_high")
    for scale in SCALES:
        sub = wide2.loc[scale]
        delta = (sub["human only"] - sub["pooled (as distributed)"]).to_numpy()
        print(f"  {scale:>7s}  human-only - pooled: {delta.mean():+7.3f}  "
              f"p={wilcoxon(delta).pvalue:.2g}  better in "
              f"{int((delta < 0).sum())}/{len(delta)}")


if __name__ == "__main__":
    main()
