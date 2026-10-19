"""Species composition of `PPBR_AZ`, and the floor it puts under the metric.

The deposition behind the benchmark holds five plasma protein binding assays,
one per species. The benchmark ships one endpoint. This module establishes which
assay each benchmark row came from, by the same value match used for hepatocyte
clearance, and writes the result so the figures in the paper are derived rather
than transcribed.

The duplication that the pooling creates also bounds the leaderboard metric from
below: where one molecule carries several contradictory labels, no prediction
satisfies all of them. That floor is computed here too, since it is a property of
the distributed data and not of any model.

Two numbers are reported for every row: the assay it matches, and whether that
match is unique. A value that appears in two species' assays for the same
molecule cannot be attributed, and 1.4% of rows are in that position.
"""

from __future__ import annotations

import json
import logging
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from tdc.benchmark_group import admet_group

from src.data.load_tdc import DATA_ROOT, canonical_smiles
from src.models.baseline import RESULTS_ROOT

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
#: The five plasma protein binding assays of deposition CHEMBL3301361. The
#: benchmark's description names neither the deposition nor the species.
PPB_ASSAYS = {
    "CHEMBL3301365": "human",
    "CHEMBL3301366": "rat",
    "CHEMBL3301367": "dog",
    "CHEMBL3301368": "mouse",
    "CHEMBL3301369": "guinea pig",
}
UPSTREAM = RESULTS_ROOT / "ppb_all_species_upstream.csv"
COMPOSITION = RESULTS_ROOT / "ppb_species_composition.csv"
FLOOR = RESULTS_ROOT / "ppb_metric_floor.csv"


def fetch_upstream() -> pd.DataFrame:
    """Every record of all five assays, with its species."""
    rows = []
    for assay_id, species in PPB_ASSAYS.items():
        offset = 0
        while True:
            query = urllib.parse.urlencode(
                {"assay_chembl_id": assay_id, "limit": 1000, "offset": offset}
            )
            with urllib.request.urlopen(
                f"{CHEMBL}/activity.json?{query}", timeout=120
            ) as handle:
                payload = json.load(handle)
            for activity in payload["activities"]:
                if activity.get("canonical_smiles") and activity.get("standard_value"):
                    rows.append({
                        "species": species,
                        "smiles": activity["canonical_smiles"],
                        "value": float(activity["standard_value"]),
                        "rel": activity.get("standard_relation"),
                    })
            if not payload["page_meta"].get("next"):
                break
            offset += 1000
        log.info("%-11s %s", species, assay_id)
    frame = pd.DataFrame(rows)
    frame["c"] = frame["smiles"].map(canonical_smiles)
    return frame


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    upstream = pd.read_csv(UPSTREAM) if UPSTREAM.exists() else fetch_upstream()
    if not UPSTREAM.exists():
        upstream.to_csv(UPSTREAM, index=False)

    group = admet_group(path=str(DATA_ROOT / "raw"))
    benchmark = group.get("PPBR_AZ")
    rows = pd.concat([benchmark["train_val"], benchmark["test"]]).copy()
    rows["c"] = [canonical_smiles(s) for s in rows["Drug"]]
    rows["key"] = list(zip(rows["c"], rows["Y"].round(4)))

    hits = {}
    for species, block in upstream.groupby("species"):
        keys = set(zip(block["c"], block["value"].round(4)))
        hits[species] = rows["key"].isin(keys)
    matched = pd.DataFrame(hits)
    unique = matched.sum(axis=1) == 1

    order = ["human", "rat", "dog", "mouse", "guinea pig"]
    composition = pd.DataFrame({
        "species": order,
        "upstream_records": [int((upstream["species"] == s).sum()) for s in order],
        "matched_rows": [int(matched[s].sum()) for s in order],
        "unambiguous_rows": [int(matched.loc[unique, s].sum()) for s in order],
    })
    composition["unambiguous_share"] = (
        composition["unambiguous_rows"] / int(unique.sum())
    ).round(4)
    composition.to_csv(COMPOSITION, index=False)

    # The floor: the best single prediction per molecule is the median of its
    # labels under absolute error, and the mean under squared error.
    test = benchmark["test"].copy()
    test["c"] = [canonical_smiles(s) for s in test["Drug"]]
    floor = pd.DataFrame([{
        "test_rows": len(test),
        "test_molecules": int(test["c"].nunique()),
        "molecules_with_conflicting_labels": int(
            (test.groupby("c")["Y"].nunique() > 1).sum()
        ),
        "max_within_molecule_spread": float(
            test.groupby("c")["Y"].agg(lambda v: v.max() - v.min()).max()
        ),
        "mae_floor": float(
            np.mean(np.abs(test["Y"] - test.groupby("c")["Y"].transform("median")))
        ),
        "rmse_floor": float(
            np.sqrt(np.mean((test["Y"] - test.groupby("c")["Y"].transform("mean")) ** 2))
        ),
    }])
    floor.to_csv(FLOOR, index=False)

    print("\n=== PPBR_AZ species composition, by value match ===")
    print(composition.to_string(index=False))
    print(
        f"\nbenchmark rows {len(rows)}, molecules {rows['c'].nunique()}; "
        f"upstream records {len(upstream)}, molecules {upstream['c'].nunique()}"
    )
    print(
        f"unambiguous {int(unique.sum())} ({100 * unique.mean():.1f}%), "
        f"non-human share of those "
        f"{100 * matched.loc[unique, ['rat', 'dog', 'mouse', 'guinea pig']].any(axis=1).mean():.1f}%"
    )
    print("\n=== Floor under the leaderboard metric ===")
    print(floor.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
