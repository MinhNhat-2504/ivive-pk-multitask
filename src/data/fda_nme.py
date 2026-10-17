"""Step 1 of the temporal external set: FDA new molecular entities since 2017.

All five TDC source datasets were compiled by 2016 (AZ in vitro data deposited
as ChEMBL CHEMBL3301361 around 2015-16; Lombardo 2016; Obach 2008), so a
compound first approved from 2017 onward cannot have been in the training data
of any model benchmarked on these leaderboards.

This module takes the openFDA drugsfda bulk file and returns the approvals that
qualify. It writes nothing to the project's data directory that is not
reproducible from the download, and the download itself is cached.

Scope, from the reconnaissance in `notes/05`: this set can validate the in vivo
heads. It cannot validate the in vitro clearance heads, because regulatory
review documents do not report intrinsic clearance in assay units at all.
"""

from __future__ import annotations

import io
import json
import logging
import zipfile
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)

BULK_URL = "https://download.open.fda.gov/drug/drugsfda/drug-drugsfda-0001-of-0001.json.zip"
CUTOFF = "2017-01-01"
NME_DESCRIPTION = "Type 1 - New Molecular Entity"

EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"


def download_bulk(destination: Path | None = None, timeout: int = 300) -> Path:
    """Fetch the drugsfda bulk file once and cache it."""
    destination = destination or EXTERNAL_ROOT / "drugsfda.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        log.info("using cached %s", destination)
        return destination

    log.info("downloading %s", BULK_URL)
    response = requests.get(BULK_URL, timeout=timeout)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        name = archive.namelist()[0]
        destination.write_bytes(archive.read(name))
    log.info("wrote %s (%.1f MB)", destination, destination.stat().st_size / 1e6)
    return destination


def extract_nmes(path: Path, cutoff: str = CUTOFF) -> pd.DataFrame:
    """Original, approved, Type 1 NME submissions on or after the cutoff.

    One row per application. `submission_type == "ORIG"` excludes supplements,
    which would otherwise re-date a long-marketed drug to a recent label change
    and silently break the temporal guarantee.
    """
    records = json.loads(path.read_text(encoding="utf-8"))["results"]
    rows = []

    for application in records:
        number = application.get("application_number", "")
        products = application.get("products") or []
        ingredients = sorted(
            {
                ingredient.get("name", "")
                for product in products
                for ingredient in (product.get("active_ingredients") or [])
            }
        )
        brands = sorted({p.get("brand_name", "") for p in products if p.get("brand_name")})

        for submission in application.get("submissions") or []:
            if submission.get("submission_type") != "ORIG":
                continue
            if submission.get("submission_status") != "AP":
                continue
            if submission.get("submission_class_code_description") != NME_DESCRIPTION:
                continue
            date = submission.get("submission_status_date", "")
            if len(date) != 8:
                continue
            iso = f"{date[:4]}-{date[4:6]}-{date[6:]}"
            if iso < cutoff:
                continue

            rows.append(
                {
                    "application_number": number,
                    "route": "NDA" if number.startswith("NDA") else "BLA" if number.startswith("BLA") else "other",
                    "approval_date": iso,
                    "year": int(date[:4]),
                    "ingredients": "; ".join(ingredients),
                    "brand_names": "; ".join(brands),
                    "n_review_docs": sum(
                        1
                        for d in (submission.get("application_docs") or [])
                        if (d.get("type") or "").lower() == "review"
                    ),
                    "has_label_doc": any(
                        (d.get("type") or "").lower() == "label"
                        for d in (submission.get("application_docs") or [])
                    ),
                }
            )

    return pd.DataFrame(rows).sort_values(["year", "application_number"]).reset_index(drop=True)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    frame = extract_nmes(download_bulk())

    EXTERNAL_ROOT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(EXTERNAL_ROOT / "nme_since_2017.csv", index=False)

    print(f"\nNME approvals on or after {CUTOFF}: {len(frame)}")
    print("\nBy year:")
    print(frame["year"].value_counts().sort_index().to_string())
    print("\nBy route:")
    print(frame["route"].value_counts().to_string())
    nda = frame[frame["route"] == "NDA"]
    print(f"\nNDA with at least one review document: {int((nda['n_review_docs'] > 0).sum())}")
    print(f"NDA with a label document: {int(nda['has_label_doc'].sum())}")
    print(f"Wrote {EXTERNAL_ROOT / 'nme_since_2017.csv'}")


if __name__ == "__main__":
    main()
