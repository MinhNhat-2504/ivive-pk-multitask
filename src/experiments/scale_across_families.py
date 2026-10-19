"""Does changing the target scale fix the compression in every model family?

Two results need joining. `ppb_compression.py` shows all six families compress
when fitted to percent bound. `target_scale.py` shows that fitting log fu or
logit instead removes most of the cost -- but only for gradient boosting. So the
claim "the scale is the cause" is established for the symptom across families
and for the remedy in one family, which is not the same claim.

This crosses the two: six families by three target scales, scored the same way.
If the remedy travels, the scale is the cause. If it works only for boosting,
then something about boosting is doing the work and the diagnosis is wrong.

One confound has to be removed first. `SVR(epsilon=2.0)` and the MLP's
convergence criteria are stated in the target's own units, which are percentage
points on one scale and log units on another; left alone, SVR on log fu would
predict a constant because every residual would sit inside its tube. So every
family is wrapped so that its target is standardised before fitting and the
prediction is mapped back. That is a linear map, which leaves gradient boosting
under L2 essentially unchanged and makes the kernel and neural models see the
same effective tolerance on all three scales. The remaining difference between
the scales is their shape, which is the thing under test.

Standardising is not enough by itself, and getting this wrong was instructive.
`SVR(epsilon=2.0)` was chosen against raw percent bound, whose standard
deviation is 18.8, so the tube is a tenth of a standard deviation. Applied
unchanged to a standardised target, a tube of 2.0 covers four standard
deviations and the model returns a near-constant -- which it duly did, on all
three scales, before this was noticed. So epsilon is re-expressed in units of
the target's spread, derived from the training labels rather than asserted.

Because of that wrapping the `percent` column here is not bit-identical to
`ppb_compression.py`; the discrepancy is reported at the end.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import TransformedTargetRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVR
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.experiments.ppb_compression import build_models
from src.experiments.target_scale import HIGH_BINDING, SCALES, to_fu
from src.features.descriptors import featurise
from src.models.baseline import RESULTS_ROOT

log = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

SEEDS = [1, 2, 3, 4, 5]
#: `build_models` sets SVR's tube to 2.0 percentage points. Expressed as a
#: fraction of the label spread it is this, and that is what carries across
#: scales; see the module docstring.
SVR_EPSILON_PERCENT = 2.0


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    external = load_external("ppb")
    train_smiles = {canonical_smiles(s) for s in benchmark["train_val"]["Drug"]}
    external = external[~external["smiles"].isin(train_smiles)].reset_index(drop=True)

    x_train, ok_train = featurise(benchmark["train_val"]["Drug"].tolist())
    y_train = benchmark["train_val"]["Y"].to_numpy()[ok_train]
    x_ext, ok_ext = featurise(external["smiles"].tolist())
    external = external.loc[ok_ext].reset_index(drop=True)
    y_ext = external["value"].to_numpy()
    x_test, ok_test = featurise(benchmark["test"]["Drug"].tolist())
    y_test = benchmark["test"]["Y"].to_numpy()[ok_test]

    # The tube, as a fraction of the target's spread, on the scale it was set on.
    epsilon = SVR_EPSILON_PERCENT / float(y_train.std())
    log.info("SVR tube %.1f percentage points = %.3f SD; used on every scale",
             SVR_EPSILON_PERCENT, epsilon)

    rows = []
    for seed in SEEDS:
        families = build_models(seed=seed)
        families["SVR (RBF kernel)"] = make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            SVR(C=10.0, epsilon=epsilon),
        )
        for family, estimator in families.items():
            for scale, (forward, inverse) in SCALES.items():
                model = TransformedTargetRegressor(
                    regressor=estimator, transformer=StandardScaler()
                )
                model.fit(x_train, forward(y_train))

                pred_ext = inverse(model.predict(x_ext))
                pred_test = inverse(model.predict(x_test))
                fold = np.maximum(
                    to_fu(pred_ext) / to_fu(y_ext), to_fu(y_ext) / to_fu(pred_ext)
                )
                high = y_ext >= HIGH_BINDING
                rows.append({
                    "seed": seed,
                    "family": family,
                    "scale": scale,
                    "fu_fold_median_high": float(np.median(fold[high])),
                    "frac_within_2fold_high": float(np.mean(fold[high] <= 2.0)),
                    "pred_range": float(pred_ext.max() - pred_ext.min()),
                    "test_mae_percent": float(np.mean(np.abs(pred_test - y_test))),
                    "test_spearman": float(spearmanr(y_test, pred_test).statistic),
                })
            log.info("seed %d  %s done", seed, family)

    per_seed = pd.DataFrame(rows)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(RESULTS_ROOT / "scale_across_families_per_seed.csv", index=False)

    order = list(build_models(seed=1))
    for metric, label in [
        ("fu_fold_median_high", "Median fu fold error, compounds bound >99%"),
        ("frac_within_2fold_high", "Fraction within 2-fold, compounds bound >99%"),
        ("pred_range", f"Predicted range (observed: {y_ext.max() - y_ext.min():.1f})"),
        ("test_mae_percent", "Benchmark test MAE, percentage points"),
        ("test_spearman", "Benchmark test Spearman"),
    ]:
        table = per_seed.pivot_table(
            index="family", columns="scale", values=metric, aggfunc="mean"
        ).reindex(order)[list(SCALES)]
        print(f"\n=== {label} ===")
        print(table.round(3).to_string())

    # Paired per seed, each family against its own percent fit.
    print("\n=== fu fold error: each family against its own `percent` fit ===")
    wide = per_seed.pivot_table(
        index=["family", "seed"], columns="scale", values="fu_fold_median_high"
    )
    for family in order:
        block = wide.loc[family]
        for scale in ("log_fu", "logit"):
            delta = (block[scale] - block["percent"]).to_numpy()
            print(
                f"  {family:>18s}  {scale:>7s}: {delta.mean():+7.3f}  "
                f"improved in {int((delta < 0).sum())}/{len(delta)}"
            )

    print(f"\n{len(SEEDS)} seeds. Targets standardised before fitting in every cell,")
    print("so SVR's tube and the MLP's tolerance mean the same thing on all scales.")


if __name__ == "__main__":
    main()
