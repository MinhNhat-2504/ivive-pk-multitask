"""The AstraZeneca assays as deposited, with species and censoring intact.

TDC's `Clearance_Hepatocyte_AZ` merges two assays from this deposition -- human
hepatocytes (408 compounds) and *rat* hepatocytes (837) -- and documents the
split only by experiment type. Two thirds of its records carry the rat value.
It also drops ChEMBL's `standard_relation` column, so the "<" and ">"
qualifiers that mark a detection-limit reading survive upstream and vanish in
conversion.

This module goes to the source instead, keeping both the species label and the
qualifier. It is what makes the human-only experiment possible, and it is the
only way to separate a left-censored 3 from a measured 3.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"

#: assay id -> (short key, species, what it measures)
AZ_ASSAYS = {
    "CHEMBL3301372": ("cl_hep_human", "human", "hepatocyte CLint, uL/min/10^6 cells"),
    "CHEMBL3301371": ("cl_hep_rat", "rat", "hepatocyte CLint, uL/min/10^6 cells"),
    "CHEMBL3301370": ("cl_mic_human", "human", "microsome CLint, uL/min/mg"),
    "CHEMBL3301365": ("ppb_human", "human", "plasma protein binding, % bound"),
}

#: Declared working ranges, quoted from the assay descriptions in the
#: deposition itself rather than inferred from the value distribution.
DECLARED_RANGE = {
    "cl_hep_human": (3.0, 150.0),
    "cl_hep_rat": (3.0, 150.0),
    "cl_mic_human": (3.0, 150.0),
    "ppb_human": (10.0, 99.95),
}


def canonical(smiles: str | None) -> str | None:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return Chem.MolToSmiles(mol) if mol is not None else None


def fetch_assay(assay_id: str, pause: float = 0.3) -> pd.DataFrame:
    """All activities for one assay, with the relation qualifier kept."""
    rows, offset = [], 0
    while True:
        response = requests.get(
            API,
            params={
                "assay_chembl_id": assay_id,
                "limit": 1000,
                "offset": offset,
                "only": "canonical_smiles,standard_value,standard_relation,molecule_chembl_id",
            },
            timeout=120,
        )
        response.raise_for_status()
        payload = response.json()
        rows.extend(payload["activities"])
        total = payload["page_meta"]["total_count"]
        if len(rows) >= total:
            break
        offset += 1000
        time.sleep(pause)

    frame = pd.DataFrame(rows)
    frame["smiles"] = frame["canonical_smiles"].map(canonical)
    frame["value"] = pd.to_numeric(frame["standard_value"], errors="coerce")
    return frame.dropna(subset=["smiles", "value"])[
        ["molecule_chembl_id", "smiles", "value", "standard_relation"]
    ]


def load(root: Path = EXTERNAL_ROOT) -> pd.DataFrame:
    """Fetch (or reuse the cache of) all four assays as one long frame."""
    cache = root / "az_assays.csv"
    if cache.exists():
        log.info("using cached %s", cache)
        return pd.read_csv(cache)

    root.mkdir(parents=True, exist_ok=True)
    frames = []
    for assay_id, (key, species, measures) in AZ_ASSAYS.items():
        log.info("fetching %s (%s)", assay_id, key)
        frame = fetch_assay(assay_id)
        frame["assay"] = key
        frame["assay_chembl_id"] = assay_id
        frame["species"] = species
        frame["measures"] = measures
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(cache, index=False)
    log.info("wrote %s", cache)
    return combined


def censoring_flags(frame: pd.DataFrame) -> pd.DataFrame:
    """Mark censored records, by qualifier and by value.

    Both are needed and they disagree. The qualifier is authoritative where it
    exists, but 13 microsome records sit at exactly the lower limit without
    one, and the PPB assay carries no qualifiers at all despite a maximum
    exactly equal to its declared ceiling. So a value at a declared bound is
    treated as censored whether or not it is flagged, and the two detections
    are kept separate so the disagreement stays visible.
    """
    frame = frame.copy()
    bounds = frame["assay"].map(DECLARED_RANGE)
    low = bounds.map(lambda b: b[0] if isinstance(b, tuple) else float("nan"))
    high = bounds.map(lambda b: b[1] if isinstance(b, tuple) else float("nan"))

    frame["censored_by_flag"] = frame["standard_relation"].isin(["<", ">"])
    frame["censored_by_value"] = (frame["value"] <= low) | (frame["value"] >= high)
    frame["censor_side"] = pd.NA
    frame.loc[frame["value"] <= low, "censor_side"] = "left"
    frame.loc[frame["value"] >= high, "censor_side"] = "right"
    return frame


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    frame = censoring_flags(load())

    print("\n=== AstraZeneca assays as deposited ===")
    summary = (
        frame.groupby(["assay", "species"])
        .agg(
            n=("value", "size"),
            n_molecules=("smiles", "nunique"),
            flagged=("censored_by_flag", "sum"),
            at_bound=("censored_by_value", "sum"),
        )
        .reset_index()
    )
    summary["pct_at_bound"] = (100 * summary["at_bound"] / summary["n"]).round(1)
    print(summary.to_string(index=False))

    print("\n=== Where flag and value disagree ===")
    for assay, group in frame.groupby("assay"):
        flagged_only = int((group["censored_by_flag"] & ~group["censored_by_value"]).sum())
        value_only = int((~group["censored_by_flag"] & group["censored_by_value"]).sum())
        print(f"  {assay:14s} flagged but not at a bound: {flagged_only:3d}   "
              f"at a bound but unflagged: {value_only:3d}")

    human = frame[(frame["assay"] == "cl_hep_human")]["smiles"]
    rat = frame[(frame["assay"] == "cl_hep_rat")]["smiles"]
    micro = frame[(frame["assay"] == "cl_mic_human")]["smiles"]
    print("\n=== Overlap, by molecule ===")
    print(f"  human hepatocyte n molecules      : {human.nunique()}")
    print(f"  rat hepatocyte   n molecules      : {rat.nunique()}")
    print(f"  human hepatocyte & human microsome: {len(set(human) & set(micro))}")
    print(f"  human hepatocyte & rat hepatocyte : {len(set(human) & set(rat))}")


if __name__ == "__main__":
    main()
