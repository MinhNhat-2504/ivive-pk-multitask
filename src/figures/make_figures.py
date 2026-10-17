"""Figures for the paper.

Three, one per claim that a reader will not take on the strength of a table:

    1  where the baseline sits between the reproducible reference and the
       corrected noise ceiling, and how wide the seed spread is around it
    2  the label-efficiency curve, with the ablation that says what causes it
    3  the censoring that governs both the ceiling estimate and the data itself

Palette is Okabe-Ito, chosen because it survives dichromacy and greyscale
printing; it passes the lightness band, chroma floor, CVD separation,
normal-vision floor and contrast checks against a paper-white surface. Every
series is also distinguishable without colour -- by position, marker or
direct label -- so identity never rests on hue alone.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT
from src.models.baseline import RESULTS_ROOT

FIGURE_ROOT = Path(__file__).resolve().parents[2] / "paper" / "figures"

BLUE, VERMILLION, GREEN = "#0072B2", "#D55E00", "#009E73"
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d4"

mpl.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "axes.labelcolor": INK,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "figure.facecolor": "white",
    }
)


def figure_1_ceiling() -> None:
    """Baseline seed spread against the reference and the corrected ceiling."""
    per_seed = pd.read_csv(RESULTS_ROOT / "leaderboard_protocol_per_seed.csv")
    resampled = pd.read_csv(RESULTS_ROOT / "baseline_lgbm_25_per_seed.csv")

    order = ["cl_hep", "cl_mic", "vdss", "half_life"]
    labels = {
        "cl_hep": "Hepatocyte CL",
        "cl_mic": "Microsome CL",
        "vdss": "VDss",
        "half_life": "Half-life",
    }
    # MapLight+GNN as re-evaluated by the end-user audit; the top-ranked entry
    # on two of these endpoints could not be obtained, so it is not the
    # reference. Ceiling is available only for hepatocyte clearance.
    reference = {"cl_hep": 0.500, "cl_mic": 0.635, "vdss": 0.716, "half_life": 0.563}
    ceiling = {"cl_hep": 0.764}

    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    rng = np.random.default_rng(0)

    for row, key in enumerate(order):
        for offset, frame, colour, marker in (
            (0.16, per_seed, BLUE, "o"),
            (-0.16, resampled, VERMILLION, "s"),
        ):
            values = frame.loc[frame["endpoint"] == key, "value"].to_numpy()
            jitter = rng.normal(0, 0.028, values.size)
            ax.scatter(
                values, np.full(values.size, row + offset) + jitter,
                s=7, color=colour, alpha=0.55, linewidths=0, marker=marker, zorder=3,
            )
            ax.plot(
                [values.mean()], [row + offset], marker="|", markersize=13,
                color=colour, markeredgewidth=1.6, zorder=4,
            )

        ax.plot([reference[key]], [row], marker="D", markersize=5.5,
                color=INK, zorder=5, clip_on=False)
        if key in ceiling:
            ax.axvline(ceiling[key], ymin=(row + 0.1) / len(order),
                       ymax=(row + 0.9) / len(order),
                       color=GREEN, linewidth=1.8, zorder=2)
            ax.annotate(
                "noise ceiling\n(lower bound)", xy=(ceiling[key], row), xytext=(4, 0),
                textcoords="offset points", fontsize=6.5, color=GREEN,
                va="center", linespacing=1.3,
            )

    ax.set_yticks(range(len(order)), [labels[k] for k in order])
    ax.set_xlabel("Spearman correlation")
    ax.set_xlim(0.25, 0.86)
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)

    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=BLUE, markersize=4,
                   label="fixed test set (leaderboard protocol)"),
        plt.Line2D([], [], marker="s", linestyle="", color=VERMILLION, markersize=4,
                   label="resampled test set"),
        plt.Line2D([], [], marker="D", linestyle="", color=INK, markersize=4.5,
                   label="reproducible reference (MapLight+GNN)"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, 1.22),
              ncol=1, fontsize=6.5, handletextpad=0.4)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "fig1_ceiling_and_spread.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "fig1_ceiling_and_spread.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_2_label_efficiency() -> None:
    """The label-efficiency curve and the ablation that explains it."""
    summary = pd.read_csv(RESULTS_ROOT / "label_efficiency_summary.csv")
    per_seed = pd.read_csv(RESULTS_ROOT / "label_efficiency_per_seed.csv")

    counts = per_seed.groupby("fraction")["n_hep_labels"].mean()
    arms = [
        ("single_task", "single-task", VERMILLION, "s", "--"),
        ("constrained", "constrained (learned coupling)", BLUE, "o", "-"),
        ("constrained_identity", "constrained (coupling fixed)", GREEN, "^", ":"),
    ]

    fig, (left, right) = plt.subplots(1, 2, figsize=(6.6, 2.7))

    for column, label, colour, marker, style in arms:
        means = per_seed.groupby("fraction")[column].mean()
        sem = per_seed.groupby("fraction")[column].sem()
        left.plot(counts, means, marker=marker, markersize=4.5, color=colour,
                  linewidth=1.6, linestyle=style, label=label, zorder=3)
        left.fill_between(counts, means - sem, means + sem, color=colour,
                          alpha=0.14, linewidth=0, zorder=2)

    left.set_xlabel("hepatocyte clearance labels used")
    left.set_ylabel("Spearman on held-out test")
    left.grid(axis="y", zorder=0)
    left.set_axisbelow(True)
    left.legend(fontsize=6.3, loc="lower right")

    # Right panel: the paired difference, which is what the statistics are on.
    deltas = summary[summary["comparison"] == "constrained - single_task"]
    ablation = summary[summary["comparison"] == "constrained - constrained_identity"]
    width = 0.36
    positions = np.arange(len(deltas))

    for shift, frame, colour, label in (
        (-width / 2, deltas, BLUE, "vs single-task"),
        (width / 2, ablation, GREEN, "vs fixed coupling"),
    ):
        frame = frame.sort_values("fraction")
        errors = np.vstack([
            frame["mean_delta"] - frame["ci_lo"],
            frame["ci_hi"] - frame["mean_delta"],
        ])
        right.bar(positions + shift, frame["mean_delta"], width, color=colour,
                  label=label, zorder=3, linewidth=0)
        right.errorbar(positions + shift, frame["mean_delta"], yerr=errors,
                       fmt="none", ecolor=INK, elinewidth=0.9, capsize=2.5, zorder=4)

    right.axhline(0, color=INK, linewidth=0.8, zorder=2)
    right.set_xticks(positions, [f"{int(f * 100)}%" for f in sorted(deltas["fraction"])])
    right.set_xlabel("fraction of labels retained")
    right.set_ylabel("paired Δ Spearman")
    right.grid(axis="y", zorder=0)
    right.set_axisbelow(True)
    right.legend(fontsize=6.3, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "fig2_label_efficiency.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "fig2_label_efficiency.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_3_censoring() -> None:
    """Where the clearance labels actually sit."""
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5), sharey=True)

    for ax, (name, key, colour) in zip(
        axes,
        (("Clearance_Hepatocyte_AZ", "hepatocyte", BLUE),
         ("Clearance_Microsome_AZ", "microsome", VERMILLION)),
    ):
        values = ADME(name=name, path=str(DATA_ROOT / "raw")).get_data()["Y"].to_numpy()
        at_low = int((values == 3.0).sum())
        at_high = int((values == 150.0).sum())
        interior = values[(values != 3.0) & (values != 150.0)]

        bins = np.logspace(np.log10(2.8), np.log10(165), 42)
        ax.hist(interior, bins=bins, color=MUTED, alpha=0.55, linewidth=0, zorder=2)
        ax.bar([3.0], [at_low], width=0.55, color=colour, zorder=3)
        ax.bar([150.0], [at_high], width=26, color=colour, zorder=3)

        ax.set_xscale("log")
        ax.set_xlabel("intrinsic clearance (µL/min)")
        ax.set_title(
            f"{key} — {100 * (at_low + at_high) / values.size:.0f}% on a limit",
            color=INK,
        )
        for x, n, align in ((3.0, at_low, "left"), (150.0, at_high, "right")):
            ax.annotate(
                f"{'<' if align == 'left' else '>'}{x:.0f}\nn={n}",
                xy=(x, n), xytext=(5 if align == "left" else -5, 2),
                textcoords="offset points", fontsize=6.5, color=colour,
                ha=align, va="bottom", linespacing=1.3,
            )
        ax.grid(axis="y", zorder=0)
        ax.set_axisbelow(True)
        # Headroom for the limit annotations, which otherwise ride into the
        # title on the panel whose lower-limit spike is tallest.
        ax.set_ylim(0, max(at_low, at_high) * 1.26)

    axes[0].set_ylabel("molecules")
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "fig3_censoring.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "fig3_censoring.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    figure_1_ceiling()
    figure_2_label_efficiency()
    figure_3_censoring()
    for path in sorted(FIGURE_ROOT.glob("*.png")):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
