"""The calibration failure that censored training produces, and its cost in fu.

Section 4.5's result needs a figure because the table understates it. A model
trained on `PPBR_AZ` ranks the external compounds better than it ranks the
benchmark's own test set, so nothing is obviously wrong until the predictions
are plotted against truth — at which point they are visibly compressed into the
interior, with the two bounds of the assay empty.

The right panel converts that into the quantity that matters downstream. The
unbound fraction is 1 - PPB/100, so a few percentage points of compression at
the top of the binding scale is a large multiplicative error in fu, and fu is
what every well-stirred clearance calculation consumes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from lightgbm import LGBMRegressor
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.features.descriptors import featurise
from src.figures.make_figures import BLUE, FIGURE_ROOT, GRID, INK, MUTED, VERMILLION
from src.models.baseline import LGBM_PARAMS

#: The assay's declared working range, from the deposition's own description.
PPB_FLOOR, PPB_CEILING = 10.0, 99.95


def build_predictions(seed: int = 1):
    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    external = load_external("ppb")

    train_smiles = {canonical_smiles(s) for s in benchmark["train_val"]["Drug"]}
    external = external[~external["smiles"].isin(train_smiles)].reset_index(drop=True)

    x_train, ok_train = featurise(benchmark["train_val"]["Drug"].tolist())
    y_train = benchmark["train_val"]["Y"].to_numpy()[ok_train]
    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)

    model = LGBMRegressor(random_state=seed, **LGBM_PARAMS).fit(x_train, y_train)
    external["pred"] = model.predict(x_ext)
    return external


def main() -> None:
    external = build_predictions()
    truth = external["value"].to_numpy()
    pred = external["pred"].to_numpy()
    censored = external["censoring"].isin(["gt", "lt"]).to_numpy()

    fig, (left, right) = plt.subplots(1, 2, figsize=(6.6, 2.9))

    # --- left: predicted against observed -------------------------------------
    left.plot([0, 100], [0, 100], color=MUTED, linewidth=0.9, linestyle="--", zorder=2)
    left.scatter(truth[~censored], pred[~censored], s=13, color=BLUE, alpha=0.6,
                 linewidths=0, zorder=3, label="measured")
    left.scatter(truth[censored], pred[censored], s=24, facecolors="none",
                 edgecolors=VERMILLION, linewidths=1.0, zorder=4, label="at an assay bound")
    for bound in (PPB_FLOOR, PPB_CEILING):
        left.axvline(bound, color=GRID, linewidth=1.0, zorder=1)

    # Above the identity line is over-prediction, below it under-prediction.
    # The first draft put these the wrong way round.
    left.annotate("over-predicted\n(weakly bound)", xy=(8, 72), fontsize=6.3,
                  color=INK, linespacing=1.3)
    left.annotate("under-predicted\n(highly bound)", xy=(97, 34), fontsize=6.3,
                  color=INK, ha="right", linespacing=1.3)
    left.set_xlabel("observed % bound (label)")
    left.set_ylabel("predicted % bound")
    left.set_xlim(0, 103)
    left.set_ylim(0, 103)
    left.legend(fontsize=6.3, loc="lower right")
    left.grid(zorder=0)
    left.set_axisbelow(True)

    # --- right: what that costs in unbound fraction ---------------------------
    high = truth >= 99.0
    fu_true = (100.0 - truth[high]) / 100.0
    fu_pred = (100.0 - np.clip(pred[high], None, 99.99)) / 100.0

    right.scatter(fu_true, fu_pred, s=16, color=VERMILLION, alpha=0.7,
                  linewidths=0, zorder=3)
    lims = [4e-4, 3e-1]
    right.plot(lims, lims, color=MUTED, linewidth=0.9, linestyle="--", zorder=2)
    for factor, style in ((10, ":"),):
        right.plot(lims, [v * factor for v in lims], color=GRID,
                   linewidth=0.9, linestyle=style, zorder=2)

    median_ratio = float(np.median(fu_pred / fu_true))
    right.annotate(
        f"median {median_ratio:.1f}× too high\n(n = {int(high.sum())}, bound > 99%)",
        # Low and right, where the identity line leaves the panel empty. The
        # first draft put this top-left, on top of the data.
        xy=(1.3e-2, 9e-4), fontsize=6.8, color=VERMILLION, linespacing=1.4,
    )
    right.annotate("10×", xy=(lims[0] * 1.3, lims[0] * 13), fontsize=6, color=MUTED)

    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlim(*lims)
    right.set_ylim(*lims)
    right.set_aspect("equal")
    right.set_xlabel("observed unbound fraction $f_u$")
    right.set_ylabel("predicted $f_u$")
    right.grid(zorder=0)
    right.set_axisbelow(True)

    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_ROOT / f"paperA_fig3_calibration.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIGURE_ROOT / 'paperA_fig3_calibration.png'}")
    print(f"  n external = {len(truth)}, at a bound = {int(censored.sum())}, "
          f"bound > 99% = {int(high.sum())}, median fu ratio = {median_ratio:.2f}")


if __name__ == "__main__":
    main()
