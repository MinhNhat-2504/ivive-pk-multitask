"""Assemble the dataset deposit: the files, a README, and the metadata.

The deposit is the external set and the provenance queries behind the two
findings, not the whole repository — a reader who wants the code goes to the
repository, and a reader who wants the data should get it without having to
run anything.

Every file ships with a column description, because the point of the paper is
that undocumented columns cause exactly the problems it reports.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "data" / "external"
RESULTS = ROOT / "results"
BUNDLE = ROOT / "zenodo_bundle"

#: file -> (short description, the columns worth explaining)
CONTENTS = {
    "ppb_final.csv": (
        "Verified plasma protein binding for drugs first approved from 2017, "
        "disjoint from TDC and from PKSmart's external sets.",
        {
            "value": "percent bound to human plasma protein",
            "censoring": "none | gt | lt | approx, as stated by the label",
            "source_quote": "the label sentence the value was read from",
            "review": "accept or correct, from the hand review",
            "approval_date": "FDA approval date of the original NDA",
            "parent_smiles": "largest fragment, salts removed",
        },
    ),
    "half_life_final.csv": (
        "Verified terminal half-life for the same cohort, in hours.",
        {
            "value": "terminal half-life in hours",
            "is_apparent": "true where the label reports an apparent (/F) value",
            "value_high": "upper bound where the label gave a range",
        },
    ),
    "ppb_manual_review.csv": (
        "Every non-accept judgement from the hand review of protein binding, "
        "with its reason. Accepted values are not listed here; they are the "
        "complement.",
        {
            "verdict": "reject | correct | suspect",
            "corrected_value": "replacement where the sentence contained the right value",
            "reason": "why, in one line",
        },
    ),
    "half_life_manual_review.csv": (
        "As above, for half-life.", {},
    ),
    "az_assays.csv": (
        "The four AstraZeneca assays as deposited in ChEMBL, with species and "
        "the measurement qualifier that the benchmark's conversion drops. This "
        "is the file behind Sections 4.1 and 4.2 of the paper.",
        {
            "assay": "cl_hep_human | cl_hep_rat | cl_mic_human | ppb_human",
            "species": "human or rat, from assay_organism",
            "standard_relation": "'<' or '>' where the value is at a detection limit",
        },
    ),
    "external_candidates_dedup.csv": (
        "All 255 temporally clean candidates before extraction, flagged for "
        "overlap with PKSmart.",
        {"in_pksmart": "true where the compound appears in PKSmart's external sets"},
    ),
    "pk_extractions.csv": (
        "Raw extractions before hand review, including rejected ones, so the "
        "context filter can be audited.",
        {"rejected_reason": "non-empty where context disqualified the value"},
    ),
    "nme_since_2017.csv": (
        "The 471 FDA new molecular entities the set is drawn from.", {},
    ),
    "ppb_all_species_upstream.csv": (
        "All five plasma protein binding assays of ChEMBL deposition "
        "CHEMBL3301361 as fetched, one row per deposited record, with the "
        "species that the benchmark's conversion drops. This is the file "
        "behind the finding that PPBR_AZ pools five species.",
        {
            "species": "human | rat | dog | mouse | guinea pig",
            "value": "percent bound to plasma of that species",
            "rel": "standard_relation; empty throughout for this endpoint",
            "c": "canonical SMILES, as used for the value match",
        },
    ),
    "ppb_species_composition.csv": (
        "The attribution of PPBR_AZ's rows to those five assays by matching "
        "on (canonical SMILES, value). 97.5% of rows match exactly one assay.",
        {
            "upstream_records": "records in the assay as deposited",
            "matched_rows": "benchmark rows whose value appears in that assay",
            "unambiguous_rows": "of those, the rows matching no other assay",
            "unambiguous_share": "fraction of all unambiguous rows",
        },
    ),
    "ppb_metric_floor.csv": (
        "The floor the per-species duplication puts under the leaderboard "
        "metric: where one molecule carries several disagreeing labels, no "
        "single prediction satisfies them all.",
        {
            "mae_floor": "mean absolute error no model can pass, percentage points",
            "rmse_floor": "the same for root mean squared error",
            "molecules_with_conflicting_labels": "of the test partition's 343",
        },
    ),
}


def locate(name: str) -> Path:
    """Where a deposit file lives.

    Curated data and derived results sit in different directories; the deposit
    flattens them, so the directory is found rather than declared per entry.
    """
    for directory in (EXTERNAL, RESULTS):
        candidate = directory / name
        if candidate.exists():
            return candidate
    return EXTERNAL / name


def write_readme() -> None:
    lines = [
        "# Verified external pharmacokinetic data for drugs approved from 2017,",
        "# and the AstraZeneca assay provenance behind TDC's clearance endpoints",
        "",
        "Supporting data for two papers: *Two undeclared data defects in the",
        "pharmacokinetic endpoints of the TDC ADMET benchmark*, and *A plasma",
        "protein binding benchmark that cannot support the use its endpoint",
        "is for*.",
        "",
        "Code: https://github.com/MinhNhat-2504/ivive-pk-multitask",
        "",
        "## What this is",
        "",
        "Two things, kept together because the paper uses them together.",
        "",
        "**The external set.** Plasma protein binding (n = 129) and terminal",
        "half-life (n = 76) for small molecules first approved by the FDA from",
        "2017 onward — after every TDC pharmacokinetic source dataset was",
        "compiled. Values are extracted from FDA product labels, then read",
        "individually against their own source sentence by hand; half-life was",
        "read twice. Compounds appearing in TDC or in PKSmart's released",
        "external sets are removed.",
        "",
        "**The assay provenance.** The four AstraZeneca assays as deposited in",
        "ChEMBL, retaining two fields the benchmark's conversion drops: the",
        "species of the assay, and the `<` / `>` qualifier marking a value at a",
        "detection limit.",
        "",
        "**The species composition of the binding endpoint.** All five plasma",
        "protein binding assays of the same deposition, human, rat, dog, mouse",
        "and guinea pig, with the row-level attribution of `PPBR_AZ` to them.",
        "42.4% of that benchmark is not human plasma, and the duplication this",
        "creates puts a floor of 1.96 percentage points under its leaderboard",
        "metric.",
        "",
        "## A caution on use",
        "",
        "Roughly a fifth of the binding values are censored — reported as",
        "`> 99%` or similar rather than measured. The `censoring` column says",
        "which. Treating those as point measurements is the error the",
        "accompanying paper is about.",
        "",
        "Half-life values marked `is_apparent` come from oral studies and are",
        "not comparable with intravenous terminal half-life; they are supplied",
        "but excluded from the verified set.",
        "",
        "## Files",
        "",
    ]
    for name, (description, columns) in CONTENTS.items():
        path = locate(name)
        rows = len(pd.read_csv(path)) if path.exists() else 0
        lines += [f"### `{name}` — {rows} rows", "", description, ""]
        if columns:
            lines += ["| column | meaning |", "|---|---|"]
            lines += [f"| `{c}` | {m} |" for c, m in columns.items()]
            lines += [""]

    lines += [
        "## Licence",
        "",
        "CC BY 4.0. The underlying FDA label text is a US government work; the",
        "ChEMBL deposition is CC BY-SA 3.0 (ChEMBL CHEMBL3301361). This deposit",
        "covers the curation, verification and provenance annotation.",
        "",
        "## Citation",
        "",
        "`[TODO: paper citation once accepted]`",
        "",
    ]
    (BUNDLE / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if BUNDLE.exists():
        shutil.rmtree(BUNDLE)
    BUNDLE.mkdir(parents=True)

    missing = []
    for name in CONTENTS:
        source = locate(name)
        if source.exists():
            shutil.copy2(source, BUNDLE / name)
        else:
            missing.append(name)

    write_readme()

    print(f"bundle at {BUNDLE}")
    for path in sorted(BUNDLE.iterdir()):
        print(f"  {path.name:34s} {path.stat().st_size / 1024:8.1f} KB")
    if missing:
        print("\nMISSING, regenerate before depositing:")
        for name in missing:
            print(f"  ! {name}")


if __name__ == "__main__":
    main()
