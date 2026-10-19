"""What an error in unbound fraction costs in predicted hepatic clearance.

Paper A establishes that a predicted binding value gives an unbound fraction
wrong by a median factor of 8.3 for highly bound compounds. That number is only
interesting if it survives the well-stirred model, and it does not survive it
uniformly: hepatic clearance is near-linear in fu while fu*CLint is small
against hepatic blood flow, and saturates at Q_h once it is not. So the cost of
an fu error depends on the compound's intrinsic clearance, and the dependence
is worth reporting rather than asserting a single figure.

Two analyses.

    parametric  across the assay's working range of CLint, what fold error in
                CL_h does a given fold error in fu produce?
    empirical   for the external compounds, combining their measured and
                predicted binding with a CLint predicted from the clean human
                microsome assay, what is the distribution of CL_h error?

The second needs a clearance prediction for compounds that have no measured
one, which is itself the situation the IVIVE chain is built for.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor

from src.data.chembl_az import DECLARED_RANGE, load
from src.experiments.external_validation import load_external
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT
from src.models.ivive import (
    HEPATIC_BLOOD_FLOW,
    MICROSOME_SCALE,
    unbound_fraction,
    well_stirred_clearance,
)

log = logging.getLogger(__name__)

#: Fold errors in fu to report. 8.3 is the median measured in paper A.
FOLD_ERRORS = [2.0, 4.0, 8.3, 16.0]


def parametric() -> pd.DataFrame:
    """How an fu error propagates, as a function of intrinsic clearance."""
    clint_in_vitro = np.array([3, 10, 30, 50, 100, 150], dtype=float)
    clint_scaled = clint_in_vitro * MICROSOME_SCALE

    rows = []
    for fu_true in (0.001, 0.005, 0.01, 0.05):
        for fold in FOLD_ERRORS:
            fu_pred = min(fu_true * fold, 1.0)
            true = well_stirred_clearance(clint_scaled, np.full_like(clint_scaled, fu_true))
            predicted = well_stirred_clearance(clint_scaled, np.full_like(clint_scaled, fu_pred))
            ratio = predicted / true
            rows.append(
                {
                    "fu_true": fu_true,
                    "fu_fold_error": fold,
                    **{
                        f"CLh_fold_error_at_CLint_{int(c)}": float(r)
                        for c, r in zip(clint_in_vitro, ratio)
                    },
                }
            )
    return pd.DataFrame(rows)


def empirical() -> pd.DataFrame:
    """The same propagation on the external compounds, end to end."""
    az = load()
    micro = (
        az[az["assay"] == "cl_mic_human"]
        .groupby("smiles", as_index=False)["value"]
        .mean()
    )
    x_micro, ok_micro = featurise(micro["smiles"].tolist())
    y_micro = np.log10(np.clip(micro["value"].to_numpy()[ok_micro], 1e-6, None))
    clearance_model = LGBMRegressor(random_state=1, **LGBM_PARAMS).fit(x_micro, y_micro)

    predictions = pd.read_csv(RESULTS_ROOT / "ppb_compression_predictions.csv")
    external = load_external("ppb")
    frame = predictions.merge(external[["smiles", "value"]], on="smiles", how="inner")

    x_ext, ok_ext = featurise(frame["smiles"].tolist())
    frame = frame.loc[ok_ext].reset_index(drop=True)
    clint = 10.0 ** clearance_model.predict(x_ext)
    clint_scaled = clint * MICROSOME_SCALE

    model_columns = [
        c for c in predictions.columns
        if c not in {"observed", "smiles", "censoring"}
    ]
    rows = []
    for column in model_columns:
        fu_true = unbound_fraction(frame["observed"].to_numpy())
        fu_pred = unbound_fraction(frame[column].to_numpy())
        true = well_stirred_clearance(clint_scaled, fu_true)
        predicted = well_stirred_clearance(clint_scaled, fu_pred)
        ratio = predicted / true

        high = frame["observed"].to_numpy() >= 99.0
        rows.append(
            {
                "model": column,
                "CLh_fold_error_median_all": float(np.median(np.maximum(ratio, 1 / ratio))),
                "CLh_fold_error_median_above99": float(
                    np.median(np.maximum(ratio[high], 1 / ratio[high]))
                ),
                "CLh_fold_error_p90_above99": float(
                    np.percentile(np.maximum(ratio[high], 1 / ratio[high]), 90)
                ),
                "frac_above_2fold": float(np.mean(np.maximum(ratio, 1 / ratio) > 2.0)),
                "n": int(len(ratio)),
                "n_above99": int(high.sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    table = parametric()
    table.to_csv(RESULTS_ROOT / "fu_propagation_parametric.csv", index=False)
    print("\n=== An fu error of 8.3-fold, propagated through the well-stirred model ===")
    subset = table[table["fu_fold_error"] == 8.3]
    columns = [c for c in table.columns if c.startswith("CLh_fold")]
    print(subset[["fu_true", *columns]].round(2).to_string(index=False))
    print(
        f"\nColumns are CL_h fold error at the stated in vitro CLint "
        f"(uL/min/mg). Hepatic blood flow is {HEPATIC_BLOOD_FLOW} mL/min/kg; "
        "the error saturates as fu*CLint approaches it."
    )

    empirical_table = empirical()
    empirical_table.to_csv(RESULTS_ROOT / "fu_propagation_empirical.csv", index=False)
    print("\n=== On the external compounds, end to end ===")
    print(empirical_table.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
