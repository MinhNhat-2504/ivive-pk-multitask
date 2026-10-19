"""Add the plasma protein binding species files to the dataset deposit.

The deposit was assembled for the first paper and holds the external set and
the four AstraZeneca assays. The second paper rests on three further files that
are not in it: the five species' binding assays as fetched, the row-level
attribution of the benchmark to those assays, and the floor the resulting
duplication puts under the leaderboard metric. A paper that cites a DOI for its
data should not need a reader to re-run anything to get it.

Those three are derived results rather than curated data, so they live in
`results/` while everything already deposited lives in `data/external/`. Giving
every entry an explicit source directory means editing each tuple of a Python
literal by string surgery, which broke the file on the first attempt. The
lookup searches both directories by name instead, and the entries stay as they
are.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUNDLER = ROOT / "src" / "data" / "make_zenodo_bundle.py"

NEW_ENTRIES = '''    "ppb_all_species_upstream.csv": (
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
}'''

LOOKUP = '''def locate(name: str) -> Path:
    """Where a deposit file lives.

    Curated data and derived results sit in different directories; the deposit
    flattens them, so the directory is found rather than declared per entry.
    """
    for directory in (EXTERNAL, RESULTS):
        candidate = directory / name
        if candidate.exists():
            return candidate
    return EXTERNAL / name


'''

README_LINES = '''        "**The species composition of the binding endpoint.** All five plasma",
        "protein binding assays of the same deposition, human, rat, dog, mouse",
        "and guinea pig, with the row-level attribution of `PPBR_AZ` to them.",
        "42.4% of that benchmark is not human plasma, and the duplication this",
        "creates puts a floor of 1.96 percentage points under its leaderboard",
        "metric.",
        "",
        "## A caution on use",'''

EDITS = [
    ('BUNDLE = ROOT / "zenodo_bundle"',
     'RESULTS = ROOT / "results"\nBUNDLE = ROOT / "zenodo_bundle"'),
    ('''    "nme_since_2017.csv": (
        "The 471 FDA new molecular entities the set is drawn from.", {},
    ),
}''',
     '''    "nme_since_2017.csv": (
        "The 471 FDA new molecular entities the set is drawn from.", {},
    ),
''' + NEW_ENTRIES),
    ("def write_readme() -> None:", LOOKUP + "def write_readme() -> None:"),
    ("        path = EXTERNAL / name", "        path = locate(name)"),
    ("        source = EXTERNAL / name", "        source = locate(name)"),
    ('        "## A caution on use",', README_LINES),
    ('''        "Supporting data for: *Two undeclared data defects in the pharmacokinetic",
        "endpoints of the TDC ADMET benchmark*.",''',
     '''        "Supporting data for two papers: *Two undeclared data defects in the",
        "pharmacokinetic endpoints of the TDC ADMET benchmark*, and *A plasma",
        "protein binding benchmark that cannot support the use its endpoint",
        "is for*.",'''),
]


def main() -> None:
    text = BUNDLER.read_text(encoding="utf-8")
    if "ppb_species_composition.csv" in text:
        print("already extended")
        return
    for old, new in EDITS:
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique or missing:\n  {old[:80]}")
        text = text.replace(old, new, 1)
    BUNDLER.write_text(text, encoding="utf-8")
    print("bundle assembler extended with the three species files")


if __name__ == "__main__":
    main()
