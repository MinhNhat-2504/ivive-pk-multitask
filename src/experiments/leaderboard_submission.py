"""What the target scale is worth under the leaderboard's own protocol.

Every comparison so far trained on the whole `train_val` partition and varied
only the model seed. That is not how the benchmark is scored. The official
protocol repartitions train and validation on each of five seeds, holds the test
partition fixed, and reports the mean and standard deviation of the metric --
mean absolute error on percent bound for this endpoint -- through the group's
own `evaluate`.

The distinction matters for the claim being made. Saying that a target
transformation would move a model on the public leaderboard requires producing
the number the leaderboard would compute, not a number computed a similar way.

Three scales are run under that protocol. Five seeds give the figure that is
comparable to a published entry; twenty-five give a tighter estimate of where
the scales actually sit, since five seeds of this benchmark carry enough spread
to hide a difference of a few tenths.

No tuning is done and no validation-based early stopping is used, so this is a
floor on what the transformation is worth rather than an attempt at a strong
entry.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT
from src.experiments.target_scale import SCALES
from src.features.descriptors import featurise
from src.models.baseline import LGBM_PARAMS, RESULTS_ROOT

log = logging.getLogger(__name__)

BENCHMARK = "PPBR_AZ"
#: The protocol's own seed set. The leaderboard reports mean and sd over these.
OFFICIAL_SEEDS = [1, 2, 3, 4, 5]
EXTENDED_SEEDS = list(range(1, 26))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get(BENCHMARK)
    test = benchmark["test"]
    x_test, ok_test = featurise(test["Drug"].tolist())

    rows = []
    for seed in EXTENDED_SEEDS:
        train, _valid = group.get_train_valid_split(
            benchmark=BENCHMARK, split_type="default", seed=seed
        )
        x_train, ok_train = featurise(train["Drug"].tolist())
        y_train = train["Y"].to_numpy()[ok_train]

        for scale, (forward, inverse) in SCALES.items():
            model = LGBMRegressor(random_state=seed, **LGBM_PARAMS)
            model.fit(x_train, forward(y_train))
            # `evaluate` compares against `test.Y` positionally, so the array
            # must be full length and in the file's order. Any molecule RDKit
            # could not featurise is filled with the training median rather
            # than dropped, which keeps the alignment and is the prediction a
            # submission would have to make anyway.
            prediction = np.full(len(test), float(np.median(y_train)))
            prediction[ok_test] = inverse(model.predict(x_test))

            result = group.evaluate(
                {BENCHMARK.lower(): prediction}, save_dict=False
            )
            metric, value = next(iter(result[BENCHMARK.lower()].items()))
            rows.append({"seed": seed, "scale": scale, "metric": metric, "value": value})
        log.info("seed %2d done", seed)

    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(RESULTS_ROOT / "leaderboard_submission_per_seed.csv", index=False)
    metric = per_seed["metric"].iloc[0]

    print(f"\n=== {BENCHMARK} under the official protocol, metric = {metric} ===")
    for label, seeds in (("official 5 seeds", OFFICIAL_SEEDS),
                         ("25 seeds", EXTENDED_SEEDS)):
        block = per_seed[per_seed["seed"].isin(seeds)]
        print(f"\n{label}")
        for scale in SCALES:
            values = block.loc[block["scale"] == scale, "value"].to_numpy()
            print(f"  {scale:>8s}  {values.mean():.3f} +/- {values.std(ddof=1):.3f}")

    block = per_seed[per_seed["seed"].isin(EXTENDED_SEEDS)]
    wide = block.pivot(index="seed", columns="scale", values="value")
    print("\npaired against `percent`, 25 seeds")
    for scale in ("log_fu", "logit"):
        delta = (wide[scale] - wide["percent"]).to_numpy()
        print(f"  {scale:>8s}: {delta.mean():+.3f} "
              f"[{np.percentile(delta, 2.5):+.3f}, {np.percentile(delta, 97.5):+.3f}]  "
              f"better in {int((delta < 0).sum())}/{len(delta)}")


if __name__ == "__main__":
    main()
