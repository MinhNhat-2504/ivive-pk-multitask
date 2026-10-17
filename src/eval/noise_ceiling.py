"""Estimate how well any model could possibly score, given assay noise.

TDC ships several measurements of the same molecule within a single endpoint.
Once SMILES are canonicalised those collapse into replicate groups, and the
spread inside a group is a direct read on the repeatability of the assay.

From that spread we get a ceiling. If the label a model is scored against is
itself a noisy draw, then even a model that recovers the true value exactly
cannot correlate perfectly with the label. Two estimates are reported:

* the classical correction-for-attenuation bound, sqrt(ICC), which applies to
  Pearson correlation; and
* a simulation for Spearman, which is the metric the leaderboard actually uses
  and for which no closed form applies.

Caveat handled explicitly below: a duplicated *record* is not a replicate
*measurement*. Groups whose values are bit-identical are counted and reported
separately, because including them would inflate the apparent repeatability.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from tdc.single_pred import ADME

from src.data.load_tdc import DATA_ROOT, PK_ENDPOINTS, canonical_smiles
from src.models.baseline import ENDPOINT_SPEC, _forward

log = logging.getLogger(__name__)

N_SIMULATIONS = 200
#: Fewer simulations inside each bootstrap resample: the resampling spread, not
#: the simulation spread, dominates the width of the interval.
BOOT_SIMULATIONS = 20
N_BOOTSTRAP = 500
RNG_SEED = 0


def replicate_groups(tdc_name: str, key: str) -> pd.DataFrame:
    """Canonicalise one endpoint's raw records without collapsing duplicates."""
    frame = ADME(name=tdc_name, path=str(DATA_ROOT / "raw")).get_data()
    frame = frame.rename(columns={"Drug": "smiles", "Y": "value"})[["smiles", "value"]]
    frame["smiles"] = frame["smiles"].map(canonical_smiles)
    frame = frame.dropna(subset=["smiles"])
    # Work in the space the model is fitted in, so the ceiling is expressed in
    # the same units as the reported metric.
    frame["value"] = _forward(frame["value"].to_numpy(), ENDPOINT_SPEC[key]["log"])
    return frame


def describe_replicates(frame: pd.DataFrame) -> dict:
    """Split replicate groups into genuine re-measurements and record copies."""
    grouped = frame.groupby("smiles")["value"]
    sizes = grouped.size()
    multi = sizes[sizes > 1].index

    identical, genuine = [], []
    for smiles in multi:
        values = frame.loc[frame["smiles"] == smiles, "value"].to_numpy()
        (identical if np.ptp(values) == 0 else genuine).append(smiles)

    return {
        "n_records": int(len(frame)),
        "n_molecules": int(sizes.size),
        "n_replicated": int(len(multi)),
        "n_identical_value_groups": int(len(identical)),
        "n_genuine_replicate_groups": int(len(genuine)),
        "genuine_smiles": genuine,
    }


def variance_components(frame: pd.DataFrame, genuine: list[str]) -> dict:
    """One-way random-effects decomposition over the genuine replicate groups.

    within  : pooled variance of measurements about their own molecule's mean
    between : variance of the molecule means, with the within-component removed
              so it estimates true biological spread rather than spread plus noise
    """
    subset = frame[frame["smiles"].isin(genuine)]
    grouped = subset.groupby("smiles")["value"]

    # Pooled within-group variance, ddof=1 per group.
    ss_within = float(((grouped.transform("mean") - subset["value"]) ** 2).sum())
    df_within = int(len(subset) - grouped.ngroups)
    if df_within <= 0:
        return {}
    var_within = ss_within / df_within

    # Between-molecule variance is estimated over *all* molecules in the
    # endpoint, not only the replicated ones: the replicated subset is small and
    # need not span the same chemical range as the endpoint as a whole.
    all_means = frame.groupby("smiles")["value"].mean()
    var_observed_between = float(all_means.var(ddof=1))
    var_between = max(var_observed_between - var_within, 0.0)

    icc = var_between / (var_between + var_within) if var_between + var_within else np.nan
    return {
        "n_groups_used": int(grouped.ngroups),
        "n_measurements_used": int(len(subset)),
        "sd_within": float(np.sqrt(var_within)),
        "sd_between": float(np.sqrt(var_between)),
        "icc": float(icc),
        "ceiling_pearson": float(np.sqrt(icc)) if np.isfinite(icc) else np.nan,
    }


def simulate_spearman_ceiling(
    frame: pd.DataFrame,
    sd_within: float,
    sd_between: float,
    n_sim: int = N_SIMULATIONS,
) -> float:
    """Highest Spearman an exactly-correct model could score against noisy labels.

    Spearman has no closed-form attenuation correction, so it is simulated. A
    model that recovers the true value exactly is correlated against a single
    noisy observation of it; that correlation bounds every model on the endpoint.

    The observed molecule means cannot serve as truth directly. Most molecules
    are measured once, so their "mean" is one noisy draw, and the spread of
    those means is inflated by exactly the noise we are trying to measure --
    correlating that against a fresh noisy draw overstates the ceiling badly.
    The means are therefore shrunk about the grand mean until their variance
    matches the de-noised between-molecule variance, which keeps the real shape
    of the distribution (log clearance is far from Gaussian) while removing the
    inflation.
    """
    means = frame.groupby("smiles")["value"].mean().to_numpy()
    observed_var = float(means.var(ddof=1))
    if observed_var <= 0:
        return float("nan")

    shrink = np.sqrt(sd_between**2 / observed_var)
    truth = means.mean() + (means - means.mean()) * shrink

    rng = np.random.default_rng(RNG_SEED)
    scores = [
        spearmanr(truth, truth + rng.normal(0.0, sd_within, truth.size)).statistic
        for _ in range(n_sim)
    ]
    return float(np.mean(scores))


def bootstrap_ceiling(
    frame: pd.DataFrame, genuine: list[str], n_boot: int = N_BOOTSTRAP
) -> dict:
    """Percentile confidence interval for the ICC and the Spearman ceiling.

    Replicate *groups* are resampled, not individual measurements: the group is
    the unit that carries the information about repeatability, and resampling
    measurements within a fixed set of groups would understate the uncertainty
    that comes from having only a few groups to begin with.
    """
    rng = np.random.default_rng(RNG_SEED)
    groups = np.asarray(genuine)

    # Per-group sum of squares and degrees of freedom, computed once. A
    # bootstrap draw is then just a sum over the drawn groups, which keeps
    # thousands of resamples cheap.
    per_group = {
        smiles: (float(((values - values.mean()) ** 2).sum()), values.size - 1)
        for smiles, values in frame[frame["smiles"].isin(groups)]
        .groupby("smiles")["value"]
        .apply(np.asarray)
        .items()
    }
    # The between-molecule term is estimated over the whole endpoint, so it is
    # held fixed across resamples; only the within term is bootstrapped.
    var_observed_between = float(frame.groupby("smiles")["value"].mean().var(ddof=1))

    iccs, ceilings = [], []
    for _ in range(n_boot):
        drawn = rng.choice(groups, size=groups.size, replace=True)
        ss = sum(per_group[s][0] for s in drawn)
        df = sum(per_group[s][1] for s in drawn)
        if df <= 0:
            continue
        var_within = ss / df
        var_between = max(var_observed_between - var_within, 0.0)
        total = var_between + var_within
        if total <= 0:
            continue
        iccs.append(var_between / total)
        ceilings.append(
            simulate_spearman_ceiling(
                frame, np.sqrt(var_within), np.sqrt(var_between), n_sim=BOOT_SIMULATIONS
            )
        )

    if not iccs:
        return {}
    return {
        "icc_lo": float(np.percentile(iccs, 2.5)),
        "icc_hi": float(np.percentile(iccs, 97.5)),
        "ceiling_spearman_lo": float(np.percentile(ceilings, 2.5)),
        "ceiling_spearman_hi": float(np.percentile(ceilings, 97.5)),
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    rows = []

    for tdc_name, key in PK_ENDPOINTS.items():
        frame = replicate_groups(tdc_name, key)
        counts = describe_replicates(frame)
        genuine = counts.pop("genuine_smiles")

        row = {"endpoint": key, **counts}
        if counts["n_genuine_replicate_groups"] >= 10:
            row.update(variance_components(frame, genuine))
            row["ceiling_spearman"] = simulate_spearman_ceiling(
                frame, row["sd_within"], row["sd_between"]
            )
            log.info("%s: bootstrapping %d resamples", key, N_BOOTSTRAP)
            row.update(bootstrap_ceiling(frame, genuine))
        else:
            log.warning(
                "%s: only %d genuine replicate groups, too few to estimate a ceiling",
                key,
                counts["n_genuine_replicate_groups"],
            )
        rows.append(row)

    report = pd.DataFrame(rows)

    print("\n=== Replicate structure of the raw data ===")
    print(
        report[
            [
                "endpoint",
                "n_records",
                "n_molecules",
                "n_replicated",
                "n_identical_value_groups",
                "n_genuine_replicate_groups",
            ]
        ].to_string(index=False)
    )

    if "icc" in report.columns:
        print("\n=== Noise ceiling (modelling scale) ===")
        cols = [
            "endpoint",
            "n_groups_used",
            "sd_within",
            "sd_between",
            "icc",
            "icc_lo",
            "icc_hi",
            "ceiling_pearson",
            "ceiling_spearman",
            "ceiling_spearman_lo",
            "ceiling_spearman_hi",
        ]
        print(report[[c for c in cols if c in report.columns]].to_string(index=False))

    out = DATA_ROOT.parent / "results" / "noise_ceiling.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
