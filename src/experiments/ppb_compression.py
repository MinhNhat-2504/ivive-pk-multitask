"""Is the compression a property of our model, or of the problem?

Paper A reports that a model trained on `PPBR_AZ` compresses its predictions
toward the interior, over-predicting weakly bound compounds and under-predicting
highly bound ones, and that this costs a median factor of 8.3 in unbound
fraction for the compounds bound above 99%. That result came from one
gradient-boosted model, and a reviewer is right to ask whether a different
model would behave differently.

So the same training data and the same features are given to six model families
with different inductive biases -- a linear model, a kernel model, an
instance-based model, two tree ensembles and a neural network -- and each is
measured the same way on the same external compounds. If all of them compress,
the conclusion is about the task and its labels rather than about any one
estimator.

The compression is quantified by the span of the predictions against the span
of the truth, by the bias in the lowest and highest binding bins, and by the
median fold error in unbound fraction, which is what propagates into clearance.
A regression slope of observed on predicted is also reported, but it is the
weakest of the four: it is attenuated by error in the predictor, so a model
that is merely inaccurate can show a slope below 1 while still spanning less
than the truth. The bin biases do not have that failure mode.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.experiments.external_validation import load_external
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=UserWarning)

#: Six families with genuinely different inductive biases. The point is the
#: spread of assumptions, not a search for the best model -- none of these is
#: tuned, because tuning each one would make the comparison about tuning.
def build_models(seed: int) -> dict:
    scaled = lambda estimator: make_pipeline(
        SimpleImputer(strategy="median"), StandardScaler(), estimator
    )
    return {
        "ridge (linear)": scaled(RidgeCV(alphas=np.logspace(-2, 4, 25))),
        "SVR (RBF kernel)": scaled(SVR(C=10.0, epsilon=2.0)),
        "k-NN (k=5)": scaled(KNeighborsRegressor(n_neighbors=5)),
        "random forest": make_pipeline(
            SimpleImputer(strategy="median"),
            RandomForestRegressor(n_estimators=400, random_state=seed, n_jobs=-1),
        ),
        "gradient boosting": LGBMRegressor(random_state=seed, **LGBM_PARAMS),
        "neural network": scaled(
            MLPRegressor(hidden_layer_sizes=(256, 64), max_iter=600, random_state=seed)
        ),
    }


def compression_metrics(truth: np.ndarray, pred: np.ndarray) -> dict:
    """How far the predictions are pulled toward the middle, three ways."""
    # Regressing observed on predicted: a calibrated model gives slope 1.
    # Below 1 means the predictions span a wider range than the truth;
    # above 1 means they are compressed relative to it.
    slope = float(np.polyfit(pred, truth, 1)[0])

    high = truth >= 99.0
    low = truth < 50.0
    fu_true = (100.0 - truth[high]) / 100.0
    fu_pred = (100.0 - np.clip(pred[high], None, 99.99)) / 100.0

    return {
        "spearman": float(spearmanr(truth, pred).statistic),
        "pred_range": float(pred.max() - pred.min()),
        "truth_range": float(truth.max() - truth.min()),
        "slope_truth_on_pred": slope,
        "bias_above_99": float(np.mean(pred[high] - truth[high])),
        "bias_below_50": float(np.mean(pred[low] - truth[low])),
        "fu_fold_error_median": float(np.median(fu_pred / fu_true)),
        "n_above_99": int(high.sum()),
    }


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

    log.info("train %d, external %d (%d bound above 99%%)",
             len(y_train), len(y_ext), int((y_ext >= 99).sum()))

    rows, predictions = [], {"observed": y_ext}
    for name, model in build_models(seed=1).items():
        model.fit(x_train, y_train)
        pred = np.clip(model.predict(x_ext), 0.0, 100.0)
        predictions[name] = pred
        rows.append({"model": name, **compression_metrics(y_ext, pred)})
        log.info("%-20s slope %.3f  bias>99%% %+.2f  fu %.1fx",
                 name, rows[-1]["slope_truth_on_pred"],
                 rows[-1]["bias_above_99"], rows[-1]["fu_fold_error_median"])

    summary = pd.DataFrame(rows)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(RESULTS_ROOT / "ppb_compression_summary.csv", index=False)
    pd.DataFrame(predictions).assign(
        smiles=external["smiles"], censoring=external["censoring"]
    ).to_csv(RESULTS_ROOT / "ppb_compression_predictions.csv", index=False)

    print("\n=== Compression across model families, same data and features ===")
    print(
        summary[["model", "spearman", "slope_truth_on_pred", "pred_range",
                 "bias_below_50", "bias_above_99", "fu_fold_error_median"]]
        .round(3).to_string(index=False)
    )
    print(f"\nobserved range: {y_ext.max() - y_ext.min():.1f} percentage points")
    print("slope of observed on predicted: 1.0 is calibrated, above 1.0 means the")
    print("predictions are compressed relative to the truth.")


if __name__ == "__main__":
    main()
