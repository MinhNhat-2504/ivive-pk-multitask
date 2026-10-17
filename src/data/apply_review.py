"""Apply the hand review to the extracted values and emit the verified sets.

Automatic extraction reached about 88% precision on protein binding and 91% on
half-life (`notes/11`), which is not good enough for a set whose purpose is to
confirm or refute a model. Every candidate was therefore read against its own
source sentence, and the judgements live in the two review files so the
curation is auditable rather than asserted.

Three verdicts. `reject` drops a value whose sentence does not contain one --
nearly always a flattened table, a bioequivalence row, or a food-effect
statement. `correct` replaces a value where the sentence does contain the right
one and the extractor took the wrong number; most of these are recoverable
because the true value sits a few words away. `suspect` keeps a value out of
the verified set without asserting it is wrong, for sentences too truncated or
too ambiguous to judge -- a parameter reported jointly for a drug and its
metabolite, or measured in an impaired population.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)
EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"

#: Each endpoint's candidate filter and the review file holding its judgements.
#: Half-life keeps only non-apparent values: an oral half-life under flip-flop
#: kinetics reflects absorption rather than elimination, and is not comparable
#: with the intravenous reference the training data uses.
ENDPOINTS = {
    "ppb": ("ppb_manual_review.csv", lambda f: f["parameter"] == "ppb"),
    "half_life": (
        "half_life_manual_review.csv",
        lambda f: (f["parameter"] == "half_life") & (~f["is_apparent"]),
    ),
}


def apply_one(extractions: pd.DataFrame, key: str, review_file: str, select) -> pd.DataFrame:
    review = pd.read_csv(EXTERNAL_ROOT / review_file)
    candidates = (
        extractions[select(extractions) & extractions["rejected_reason"].isna()]
        .drop_duplicates("application_number")
        .copy()
    )
    before = len(candidates)

    unmatched = set(review["ingredient"]) - set(candidates["ingredient"])
    if unmatched:
        log.warning("%s: review names not in extractions: %s", key, sorted(unmatched))

    candidates["review"] = (
        candidates["ingredient"]
        .map(dict(zip(review["ingredient"], review["verdict"])))
        .fillna("accept")
    )
    corrected = candidates["ingredient"].map(
        dict(zip(review["ingredient"], review["corrected_value"]))
    )
    applied = candidates["review"].eq("correct") & corrected.notna()
    candidates.loc[applied, "value"] = corrected[applied]

    verified = candidates[candidates["review"].isin(["accept", "correct"])].copy()
    verified.to_csv(EXTERNAL_ROOT / f"{key}_verified.csv", index=False)

    counts = candidates["review"].value_counts()
    accepted, fixed = int(counts.get("accept", 0)), int(counts.get("correct", 0))
    print(f"\n=== {key} ===")
    print(f"candidates : {before}")
    for verdict in ("accept", "correct", "suspect", "reject"):
        print(f"  {verdict:8s} : {int(counts.get(verdict, 0))}")
    print(f"VERIFIED   : {len(verified)}")
    print(
        f"precision  : {accepted / before:.1%} accepted outright, "
        f"{(accepted + fixed) / before:.1%} usable after correction"
    )
    print(verified["value"].describe().round(2).to_string())
    if verified["censoring"].notna().any():
        print("censoring :", verified["censoring"].fillna("none").value_counts().to_dict())
    return verified


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    extractions = pd.read_csv(EXTERNAL_ROOT / "pk_extractions.csv")
    for key, (review_file, select) in ENDPOINTS.items():
        apply_one(extractions, key, review_file, select)
    print(f"\nWrote verified sets to {EXTERNAL_ROOT}")


if __name__ == "__main__":
    main()
