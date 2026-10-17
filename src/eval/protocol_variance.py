"""How much apparent improvement can the evaluation protocol manufacture?

The TDC ADMET leaderboard reports mean +/- sd over five scaffold-split seeds.
Two questions follow, and both are answered here by resampling the 25-seed
baseline run rather than by argument:

1. Is a standard deviation estimated from five seeds trustworthy? With four
   degrees of freedom it is not, and the size of the error can be shown.

2. If one model is run twice under the protocol -- five seeds each time, the
   same model, the same data -- how large a difference in reported means can
   appear from seed choice alone? Any claimed improvement smaller than that
   difference is indistinguishable from luck.

Question 2 is the one that matters for reading the leaderboard, because the
margins between leading entries are routinely under 0.01.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.models.baseline import RESULTS_ROOT

N_DRAWS = 20_000
PROTOCOL_SEEDS = 5
RNG_SEED = 0


def sd_estimate_spread(values: np.ndarray, rng: np.random.Generator) -> dict:
    """Distribution of the sd a 5-seed run would report, given the 25-seed truth."""
    draws = np.array(
        [rng.choice(values, size=PROTOCOL_SEEDS, replace=False).std(ddof=1) for _ in range(N_DRAWS)]
    )
    return {
        "sd_25seed": float(values.std(ddof=1)),
        "sd5_median": float(np.median(draws)),
        "sd5_p2.5": float(np.percentile(draws, 2.5)),
        "sd5_p97.5": float(np.percentile(draws, 97.5)),
    }


def null_improvement(values: np.ndarray, rng: np.random.Generator) -> dict:
    """Apparent gap between two 5-seed runs of the SAME model.

    Disjoint seed sets are drawn, so the two runs share no split. The spread of
    the difference is the size of improvement the protocol can invent.
    """
    gaps = np.empty(N_DRAWS)
    for i in range(N_DRAWS):
        picked = rng.choice(values, size=2 * PROTOCOL_SEEDS, replace=False)
        gaps[i] = picked[:PROTOCOL_SEEDS].mean() - picked[PROTOCOL_SEEDS:].mean()
    absolute = np.abs(gaps)
    return {
        "null_gap_median_abs": float(np.median(absolute)),
        "null_gap_p95_abs": float(np.percentile(absolute, 95)),
        "null_gap_max_abs": float(absolute.max()),
    }


def main() -> None:
    per_seed = pd.read_csv(RESULTS_ROOT / "baseline_lgbm_25_per_seed.csv")
    rng = np.random.default_rng(RNG_SEED)
    rows = []

    for (endpoint, metric), group in per_seed.groupby(["endpoint", "metric"]):
        values = group["value"].to_numpy()
        if values.size < 2 * PROTOCOL_SEEDS:
            continue
        rows.append(
            {
                "endpoint": endpoint,
                "metric": metric,
                "n_seeds": values.size,
                "mean": float(values.mean()),
                **sd_estimate_spread(values, rng),
                **null_improvement(values, rng),
            }
        )

    report = pd.DataFrame(rows)

    print("\n=== Can a 5-seed run estimate its own variance? ===")
    print(
        report[["endpoint", "sd_25seed", "sd5_median", "sd5_p2.5", "sd5_p97.5"]]
        .round(4)
        .to_string(index=False)
    )

    print("\n=== Apparent improvement between two 5-seed runs of the SAME model ===")
    print(
        report[["endpoint", "metric", "null_gap_median_abs", "null_gap_p95_abs", "null_gap_max_abs"]]
        .round(4)
        .to_string(index=False)
    )

    out = RESULTS_ROOT / "protocol_variance.csv"
    report.to_csv(out, index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
