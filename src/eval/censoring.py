"""Censoring in the AZ clearance assays, and what it does to the ceiling estimate.

Two findings live here, and the second one corrects the first version of
`noise_ceiling.py`.

1. The AZ clearance assays are censored at both ends. Values pile up at the
   detection limits, so a large share of what is modelled as a continuous
   regression target is a limit flag rather than a measurement.

2. The replicate groups used to estimate the noise ceiling are a biased sample.
   TDC carries no exactly-duplicated (SMILES, value) record, which means
   concordant repeat measurements were collapsed upstream. Only discordant
   repeats survive as visible replicate groups, so the within-group variance
   they give is an overestimate and the resulting ceiling is a lower bound.

The sensitivity analysis below brackets how much of (2) is driven by (1).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.models.baseline import ENDPOINT_SPEC, RESULTS_ROOT, _forward
from src.eval.noise_ceiling import simulate_spearman_ceiling

#: Detection limits of the AstraZeneca intrinsic clearance assays, uL/min per
#: 10^6 cells (hepatocyte) or per mg protein (microsome). Values reported at a
#: limit are "<3" / ">150" recorded as a bare number.
ASSAY_LIMITS = {"cl_hep": (3.0, 150.0), "cl_mic": (3.0, 150.0)}


def load_raw(tdc_name: str) -> pd.DataFrame:
    frame = ADME(name=tdc_name, path=str(DATA_ROOT / "raw")).get_data()
    frame = frame.rename(columns={"Drug": "smiles", "Y": "value"})[["smiles", "value"]]
    frame["smiles"] = frame["smiles"].map(canonical_smiles)
    return frame.dropna(subset=["smiles"])


def censoring_profile(frame: pd.DataFrame, lo: float, hi: float) -> dict:
    """How much of the endpoint is detection-limit flags rather than measurement."""
    y = frame["value"].to_numpy()
    n_lo, n_hi = int((y == lo).sum()), int((y == hi).sum())
    return {
        "n_records": int(y.size),
        "n_at_lower_limit": n_lo,
        "pct_at_lower_limit": round(100 * n_lo / y.size, 1),
        "n_at_upper_limit": n_hi,
        "pct_at_upper_limit": round(100 * n_hi / y.size, 1),
        "pct_censored": round(100 * (n_lo + n_hi) / y.size, 1),
    }


def tie_structure(frame: pd.DataFrame, lo: float, hi: float) -> dict:
    """Share of rankable pairs destroyed by the two limit blocks.

    Spearman is a rank statistic. Molecules sharing the censored value carry no
    rank information relative to one another, so every pair inside a limit block
    is uninformative no matter how good the model is.
    """
    y = frame.groupby("smiles")["value"].mean().to_numpy()
    n = y.size
    total_pairs = n * (n - 1) // 2
    tied = 0
    for limit in (lo, hi):
        k = int((y == limit).sum())
        tied += k * (k - 1) // 2
    return {
        "n_molecules": n,
        "pct_pairs_tied_at_limits": round(100 * tied / total_pairs, 1),
    }


def ceiling_sensitivity(frame: pd.DataFrame, key: str, lo: float, hi: float) -> list[dict]:
    """Recompute the ceiling under progressively stricter replicate selections.

    `all discordant`   -- every visible replicate group (the original estimate)
    `uncensored only`  -- groups where no member sits at a detection limit, so
                          the spread is assay noise rather than a limit artefact
    """
    scaled = frame.assign(value=_forward(frame["value"].to_numpy(), ENDPOINT_SPEC[key]["log"]))
    limits_scaled = set(_forward(np.array([lo, hi]), ENDPOINT_SPEC[key]["log"]).tolist())

    sizes = scaled.groupby("smiles")["value"].size()
    replicated = sizes[sizes > 1].index

    subsets = {"all discordant": [], "uncensored only": []}
    for smiles in replicated:
        values = scaled.loc[scaled["smiles"] == smiles, "value"].to_numpy()
        subsets["all discordant"].append(smiles)
        if not any(v in limits_scaled for v in values):
            subsets["uncensored only"].append(smiles)

    var_observed_between = float(scaled.groupby("smiles")["value"].mean().var(ddof=1))
    rows = []
    for label, group_ids in subsets.items():
        subset = scaled[scaled["smiles"].isin(group_ids)]
        grouped = subset.groupby("smiles")["value"]
        df = int(len(subset) - grouped.ngroups)
        if df <= 0:
            continue
        var_within = float(((grouped.transform("mean") - subset["value"]) ** 2).sum()) / df
        var_between = max(var_observed_between - var_within, 0.0)
        icc = var_between / (var_between + var_within)
        rows.append(
            {
                "endpoint": key,
                "replicate_subset": label,
                "n_groups": grouped.ngroups,
                "sd_within": np.sqrt(var_within),
                "icc": icc,
                "ceiling_spearman": simulate_spearman_ceiling(
                    scaled, np.sqrt(var_within), np.sqrt(var_between)
                ),
            }
        )
    return rows


def main() -> None:
    profiles, ties, sensitivity = [], [], []

    for tdc_name, key in [
        ("Clearance_Hepatocyte_AZ", "cl_hep"),
        ("Clearance_Microsome_AZ", "cl_mic"),
    ]:
        lo, hi = ASSAY_LIMITS[key]
        frame = load_raw(tdc_name)
        profiles.append({"endpoint": key, **censoring_profile(frame, lo, hi)})
        ties.append({"endpoint": key, **tie_structure(frame, lo, hi)})
        sensitivity.extend(ceiling_sensitivity(frame, key, lo, hi))

    print("\n=== Censoring at the assay detection limits ===")
    print(pd.DataFrame(profiles).to_string(index=False))

    print("\n=== Rank information destroyed by the limit blocks ===")
    print(pd.DataFrame(ties).to_string(index=False))

    print("\n=== Ceiling sensitivity to which replicate groups are used ===")
    print(pd.DataFrame(sensitivity).round(4).to_string(index=False))
    print(
        "\nBoth rows are conditioned on discordance, because TDC retains no\n"
        "exactly-duplicated record. Every ceiling here is therefore a LOWER bound."
    )

    out = RESULTS_ROOT / "censoring.csv"
    pd.DataFrame(sensitivity).to_csv(out, index=False)
    pd.DataFrame(profiles).merge(pd.DataFrame(ties), on="endpoint").to_csv(
        RESULTS_ROOT / "censoring_profile.csv", index=False
    )
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
