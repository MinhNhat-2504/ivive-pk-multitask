"""Emit paper C's supporting tables from the results files.

The supporting information carries the full versions of tables the manuscript
gives in summary: every metric of the target-scale comparison rather than the
six that fit, every family-by-scale cell rather than one column, and the
negative control in full. Transcribing those by hand is how a number drifts, so
they are printed from the files and pasted, and the checker then requires the
printed values to still appear.

Run it and redirect into a scratch file; the tables are copied into
`paper/paper_C_supporting_information.md`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from src.models.baseline import RESULTS_ROOT

FAMILIES = ["ridge (linear)", "SVR (RBF kernel)", "k-NN (k=5)",
            "random forest", "gradient boosting", "neural network"]
SHORT = {"ridge (linear)": "Ridge", "SVR (RBF kernel)": "SVR", "k-NN (k=5)": "$k$-NN",
         "random forest": "Random forest", "gradient boosting": "Gradient boosting",
         "neural network": "Neural network"}
SCALES = ["percent", "log_fu", "logit"]


def to_markdown(frame: pd.DataFrame) -> str:
    """A markdown table, written out rather than pulling in `tabulate`."""
    columns = [str(c) for c in frame.columns]
    rows = [[("" if pd.isna(value) else str(value)) for value in row]
            for row in frame.itertuples(index=False)]
    lines = ["| " + " | ".join(columns) + " |",
             "|" + "|".join("---" for _ in columns) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def emit(title: str, frame: pd.DataFrame) -> None:
    print(f"\n**{title}**\n")
    print(to_markdown(frame))


def target_scale_full() -> None:
    per_seed = pd.read_csv(RESULTS_ROOT / "target_scale_per_seed.csv")
    metrics = [c for c in per_seed.columns if c not in {"seed", "scale"}]
    rows = []
    for metric in metrics:
        row = {"Metric": metric}
        for scale in SCALES:
            block = per_seed.loc[per_seed["scale"] == scale, metric]
            row[scale] = f"{block.mean():.4f} ± {block.std(ddof=1):.4f}"
        wide = per_seed.pivot(index="seed", columns="scale", values=metric)
        delta = (wide["logit"] - wide["percent"]).to_numpy()
        if np.ptp(delta) == 0:
            # `n_high` is a property of the evaluation set, identical on every
            # seed, so there is no paired difference to test.
            row["logit − percent"] = "constant across seeds"
        else:
            row["logit − percent"] = (
                f"{delta.mean():+.4f} [{np.percentile(delta, 2.5):+.4f}, "
                f"{np.percentile(delta, 97.5):+.4f}], p={wilcoxon(delta).pvalue:.1g}"
            )
        rows.append(row)
    emit("Table S1. Target scale, every metric, 25 seeds, trained on train_val.",
         pd.DataFrame(rows))


def families_full() -> None:
    per_seed = pd.read_csv(RESULTS_ROOT / "scale_across_families_per_seed.csv")
    metrics = ["fu_fold_median_high", "frac_within_2fold_high", "pred_range",
               "test_mae_percent", "test_spearman"]
    for metric in metrics:
        table = per_seed.pivot_table(index="family", columns="scale",
                                     values=metric, aggfunc="mean").reindex(FAMILIES)
        spread = per_seed.pivot_table(index="family", columns="scale",
                                      values=metric, aggfunc="std").reindex(FAMILIES)
        frame = pd.DataFrame({"Family": [SHORT[f] for f in FAMILIES]})
        for scale in SCALES:
            frame[scale] = [f"{m:.3f} ± {s:.3f}"
                            for m, s in zip(table[scale], spread[scale])]
        emit(f"Table S2 ({metric}). Six families by three scales, 5 seeds.", frame)


def leaderboard_full() -> None:
    per_seed = pd.read_csv(RESULTS_ROOT / "leaderboard_submission_per_seed.csv")
    rows = []
    for scale in SCALES:
        row = {"Target": scale}
        for label, seeds in (("Official 5 seeds", [1, 2, 3, 4, 5]),
                             ("25 seeds", list(range(1, 26)))):
            block = per_seed[(per_seed["scale"] == scale)
                             & (per_seed["seed"].isin(seeds))]["value"]
            row[label] = f"{block.mean():.3f} ± {block.std(ddof=1):.3f}"
        rows.append(row)
    emit("Table S3. PPBR_AZ under the group's own protocol, MAE in percentage points.",
         pd.DataFrame(rows))


def negative_control() -> None:
    per_seed = pd.read_csv(RESULTS_ROOT / "second_endpoint_per_seed.csv")
    metrics = [c for c in per_seed.columns
               if c not in {"seed", "scale", "n_interior"}]
    order = ["raw", "log10", "sqrt", "reciprocal"]
    rows = []
    for metric in metrics:
        row = {"Metric": metric}
        for scale in order:
            block = per_seed.loc[per_seed["scale"] == scale, metric]
            row[scale] = f"{block.mean():.4f} ± {block.std(ddof=1):.4f}"
        rows.append(row)
    emit("Table S4. Clearance_Microsome_AZ, target scale, 25 seeds. "
         "`log10` is the field's convention and the consumed scale.",
         pd.DataFrame(rows))


def compression_full() -> None:
    summary = pd.read_csv(RESULTS_ROOT / "ppb_compression_summary.csv")
    summary = summary.set_index("model").reindex(FAMILIES).reset_index()
    summary["model"] = [SHORT[f] for f in FAMILIES]
    columns = ["model", "spearman", "pred_range", "bias_below_50", "bias_above_99",
               "fu_overprediction_median", "fu_fold_error_median", "n_overpredicted",
               "n_at_assay_ceiling"]
    emit("Table S5. Compression on the external set under a percent-bound target, "
         "seed 1, untransformed targets.",
         summary[columns].round(3))


def species_and_floor() -> None:
    emit("Table S6. Species composition of PPBR_AZ by value match.",
         pd.read_csv(RESULTS_ROOT / "ppb_species_composition.csv"))
    emit("Table S7. The floor under the leaderboard metric.",
         pd.read_csv(RESULTS_ROOT / "ppb_metric_floor.csv").round(4).T.reset_index()
         .rename(columns={"index": "Quantity", 0: "Value"}))


def main() -> None:
    compression_full()
    target_scale_full()
    families_full()
    leaderboard_full()
    negative_control()
    species_and_floor()


if __name__ == "__main__":
    main()
