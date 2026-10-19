"""Figures for paper C.

Three, one per claim a reader will not take from a table:

    1  the compression under a percent target, and what the logit target does
       to it -- shown twice, in percentage points where it looks modest and in
       unbound fraction where it does not
    2  why the well-stirred model does not absorb the error: the damping as a
       function of intrinsic clearance, at four binding levels
    3  that the remedy travels across six model families

Palette and rcParams follow `make_figures.py`, for the same reasons: Okabe-Ito
survives dichromacy and greyscale, and every series is separable without colour.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.experiments.ppb_compression import build_models
from src.experiments.target_scale import HIGH_BINDING, SCALES, to_fu
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.ivive import HEPATIC_BLOOD_FLOW, MICROSOME_SCALE, well_stirred_clearance

FIGURE_ROOT = Path(__file__).resolve().parents[2] / "paper" / "figures"

BLUE, VERMILLION, GREEN = "#0072B2", "#D55E00", "#009E73"
PURPLE, AMBER = "#CC79A7", "#E69F00"
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


def _external_predictions(seed: int = 1) -> pd.DataFrame:
    """Percent- and logit-target predictions for the external compounds."""
    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    external = load_external("ppb")
    train_smiles = {canonical_smiles(s) for s in benchmark["train_val"]["Drug"]}
    external = external[~external["smiles"].isin(train_smiles)].reset_index(drop=True)

    x_train, ok_train = featurise(benchmark["train_val"]["Drug"].tolist())
    y_train = benchmark["train_val"]["Y"].to_numpy()[ok_train]
    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)

    out = {"observed": external["value"].to_numpy(),
           "censoring": external["censoring"].to_numpy()}
    for scale in ("percent", "logit"):
        forward, inverse = SCALES[scale]
        model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
        model.fit(x_train, forward(y_train))
        out[scale] = inverse(model.predict(x_ext))
    return pd.DataFrame(out)


def figure_1_compression() -> None:
    """The same two models, in percentage points and in unbound fraction."""
    frame = _external_predictions()
    truth = frame["observed"].to_numpy()
    high = truth >= HIGH_BINDING

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))

    # The span of the predictions is deliberately not annotated here. On this
    # seed the two models span 71 and 74 percentage points, which invites the
    # reader to conclude they are equivalent; the whole point is that they
    # differ in *where* they compress, and that the difference is only legible
    # once the axis is the unbound fraction. The shaded band marks that region
    # and each panel carries the error it produces there.
    for ax, scale, colour, title in (
        (axes[0], "percent", VERMILLION, "Percent-bound target"),
        (axes[1], "logit", BLUE, "Logit target"),
    ):
        pred = frame[scale].to_numpy()
        fu_pred = to_fu(pred[high])
        fold = np.median(np.maximum(fu_pred / to_fu(truth[high]),
                                    to_fu(truth[high]) / fu_pred))
        ax.axvspan(HIGH_BINDING, 103, color=INK, alpha=0.06, zorder=0)
        ax.plot([0, 100], [0, 100], "--", color=MUTED, lw=0.7, zorder=1)
        ax.scatter(truth[~high], pred[~high], s=9, color=colour, alpha=0.55,
                   linewidths=0, zorder=3)
        ax.scatter(truth[high], pred[high], s=16, facecolors="none",
                   edgecolors=colour, linewidths=0.8, zorder=4)
        ax.set_xlim(0, 103)
        ax.set_ylim(0, 103)
        ax.set_xlabel("Observed binding (%)")
        ax.set_title(title, color=colour)
        ax.grid(True, lw=0.4, alpha=0.6)
        ax.annotate(
            f"bound > {HIGH_BINDING:.0f}%:\nmedian $f_u$ error {fold:.1f}$\\times$",
            xy=(4, 99), fontsize=6.5, color=colour, va="top",
        )
    axes[0].set_ylabel("Predicted binding (%)")

    # The same predictions as unbound fraction, for the compounds that matter.
    ax = axes[2]
    # Equal limits on both axes, so the identity line reads at 45 degrees and
    # the asymmetry is visible: every observed value is confined to the narrow
    # band that "bound above 99%" allows, while the predictions are not.
    lim = (4e-4, 5e-1)
    ax.plot(lim, lim, "--", color=MUTED, lw=0.7, zorder=1)
    for factor in (10.0,):
        ax.plot(lim, [v * factor for v in lim], ":", color=MUTED, lw=0.6, zorder=1)
        ax.plot(lim, [v / factor for v in lim], ":", color=MUTED, lw=0.6, zorder=1)
    for scale, colour, marker, label in (
        ("percent", VERMILLION, "o", "percent target"),
        ("logit", BLUE, "s", "logit target"),
    ):
        fu_pred = to_fu(frame[scale].to_numpy()[high])
        fold = np.median(np.maximum(fu_pred / to_fu(truth[high]),
                                    to_fu(truth[high]) / fu_pred))
        ax.scatter(to_fu(truth[high]), fu_pred, s=14, color=colour, alpha=0.7,
                   linewidths=0, marker=marker, zorder=3,
                   label=f"{label} ({fold:.1f}$\\times$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("Observed $f_u$")
    ax.set_ylabel("Predicted $f_u$")
    ax.set_title(f"Unbound fraction, bound > {HIGH_BINDING:.0f}%", color=INK)
    ax.grid(True, lw=0.4, alpha=0.6)
    ax.legend(fontsize=6.5, loc="upper left")

    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "paperC_fig1_compression.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "paperC_fig1_compression.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_2_propagation() -> None:
    """The damping, and the regime in which it is unavailable."""
    clint = np.logspace(np.log10(3.0), np.log10(150.0), 200)
    scaled = clint * MICROSOME_SCALE
    fold_error = 8.3

    fig, (left, right) = plt.subplots(1, 2, figsize=(6.4, 2.6))

    levels = [
        (0.001, "bound 99.9%", BLUE, "-"),
        (0.005, "bound 99.5%", GREEN, "--"),
        (0.010, "bound 99.0%", AMBER, "-."),
        (0.050, "bound 95.0%", VERMILLION, ":"),
    ]
    for fu, label, colour, style in levels:
        true = well_stirred_clearance(scaled, np.full_like(scaled, fu))
        pred = well_stirred_clearance(scaled, np.full_like(scaled, min(fu * fold_error, 1.0)))
        left.plot(clint, pred / true, style, color=colour, lw=1.3, label=label)
        # The quantity that decides whether damping is available at all.
        right.plot(clint, 100 * fu * scaled / HEPATIC_BLOOD_FLOW, style,
                   color=colour, lw=1.3, label=label)

    left.axhline(fold_error, color=MUTED, lw=0.7, ls="--")
    left.annotate(f"the error in $f_u$ ({fold_error}$\\times$)", xy=(30, 8.85),
                  fontsize=6.5, color=MUTED, ha="center")
    left.set_xscale("log")
    left.set_xlabel("In vitro $\\mathrm{CL_{int}}$ (µL/min/mg)")
    left.set_ylabel("Fold error in $\\mathrm{CL}_h$")
    left.set_ylim(0, 9.5)
    left.set_title("The error that survives", color=INK)
    left.grid(True, lw=0.4, alpha=0.6)
    left.legend(fontsize=6.5, loc="lower left")

    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlabel("In vitro $\\mathrm{CL_{int}}$ (µL/min/mg)")
    right.set_ylabel("$f_u\\,\\mathrm{CL_{int}}$ as % of $Q_h$")
    right.set_title("Why: how far from saturation", color=INK)
    right.axhline(100, color=MUTED, lw=0.7, ls="--")
    right.annotate("saturation ($f_u\\mathrm{CL_{int}} = Q_h$)", xy=(3.4, 108),
                   fontsize=6.5, color=MUTED)
    right.grid(True, lw=0.4, alpha=0.6)

    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "paperC_fig2_propagation.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "paperC_fig2_propagation.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_3_families() -> None:
    """Six inductive biases, three target scales, one direction of effect."""
    per_seed = pd.read_csv(RESULTS_ROOT / "scale_across_families_per_seed.csv")
    order = list(build_models(seed=1))
    short = {
        "ridge (linear)": "ridge", "SVR (RBF kernel)": "SVR", "k-NN (k=5)": "$k$-NN",
        "random forest": "forest", "gradient boosting": "boosting",
        "neural network": "net",
    }
    table = per_seed.pivot_table(
        index="family", columns="scale", values="fu_fold_median_high", aggfunc="mean"
    ).reindex(order)
    spread = per_seed.pivot_table(
        index="family", columns="scale", values="fu_fold_median_high", aggfunc="std"
    ).reindex(order)

    fig, ax = plt.subplots(figsize=(5.4, 2.6))
    width = 0.26
    positions = np.arange(len(order))
    for offset, scale, colour, label in (
        (-width, "percent", VERMILLION, "percent"),
        (0.0, "log_fu", GREEN, "$\\log f_u$"),
        (width, "logit", BLUE, "logit"),
    ):
        ax.bar(positions + offset, table[scale], width * 0.92, color=colour,
               alpha=0.85, label=label, zorder=3,
               yerr=spread[scale], error_kw={"lw": 0.7, "ecolor": MUTED})

    ax.axhline(2.0, color=INK, lw=0.8, ls="--", zorder=4)
    ax.annotate("two-fold", xy=(len(order) - 0.45, 2.05), fontsize=6.5, color=INK)
    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels([short[f] for f in order])
    ax.set_ylabel("Median $f_u$ fold error\n(bound > 99%)")
    ax.set_yticks([2, 3, 5, 10, 20])
    ax.set_yticklabels(["2", "3", "5", "10", "20"])
    ax.grid(True, axis="y", lw=0.4, alpha=0.6)
    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.18))

    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "paperC_fig3_families.png", bbox_inches="tight")
    fig.savefig(FIGURE_ROOT / "paperC_fig3_families.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    figure_1_compression()
    figure_2_propagation()
    figure_3_families()
    for path in sorted(FIGURE_ROOT.glob("paperC_*.png")):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
