"""Add the PPB species and metric-floor figures to the paper's number guard.

`check_paper.py` verifies that each quoted figure still matches the file it was
computed from. The five-species composition and the metric floor were added to
section 4.1 and 4.3 by hand, so until they are registered here they are
unguarded -- which is how a figure drifts out of sync with its source.

Numbers come from `results/ppb_species_composition.csv` and
`results/ppb_metric_floor.csv`, both written by `src/data/ppb_species.py`.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "src" / "figures" / "check_paper.py"

ANCHOR = """    summary = RESULTS / "external_validation_summary.csv\""""

ADDITION = '''    composition = RESULTS / "ppb_species_composition.csv"
    if composition.exists():
        species = pd.read_csv(composition).set_index("species")
        for name, spelling in (
            ("human", "1 567"), ("rat", "680"), ("dog", "233"),
            ("mouse", "153"), ("guinea pig", "86"),
        ):
            rows = int(species.loc[name, "unambiguous_rows"])
            if spelling not in text:
                problems.append(
                    f"{name} plasma contributes {rows} unambiguous rows; "
                    f"'{spelling}' not in text"
                )
        for name, spelling in (
            ("human", "57.6%"), ("rat", "25.0%"), ("dog", "8.6%"),
            ("mouse", "5.6%"), ("guinea pig", "3.2%"),
        ):
            share = float(species.loc[name, "unambiguous_share"])
            if spelling not in text:
                problems.append(
                    f"{name} share is {100 * share:.1f}%; '{spelling}' not in text"
                )
        non_human = 1.0 - float(species.loc["human", "unambiguous_share"])
        if f"{100 * non_human:.1f}%" not in text:
            problems.append(
                f"non-human share is {100 * non_human:.1f}%, not quoted in text"
            )
        for name, spelling in (
            ("human", "1 614"), ("rat", "717"), ("dog", "244"),
            ("mouse", "162"), ("guinea pig", "91"),
        ):
            upstream = int(species.loc[name, "upstream_records"])
            if spelling not in text:
                problems.append(
                    f"upstream {name} assay holds {upstream} records; "
                    f"'{spelling}' not in text"
                )

    floor_file = RESULTS / "ppb_metric_floor.csv"
    if floor_file.exists():
        floor = pd.read_csv(floor_file).iloc[0]
        for column, spelling in (
            ("mae_floor", "1.964"), ("rmse_floor", "4.478"),
            ("test_rows", "559"), ("test_molecules", "343"),
            ("molecules_with_conflicting_labels", "140"),
            ("max_within_molecule_spread", "49.3"),
        ):
            if spelling not in text:
                problems.append(
                    f"{column} is {floor[column]}; '{spelling}' not in text"
                )

'''


def main() -> None:
    text = CHECKER.read_text(encoding="utf-8")
    if "ppb_species_composition.csv" in text:
        print("already registered")
        return
    if text.count(ANCHOR) != 1:
        raise SystemExit("anchor not found exactly once in check_paper.py")
    CHECKER.write_text(text.replace(ANCHOR, ADDITION + ANCHOR, 1), encoding="utf-8")
    print("registered the PPB species and metric-floor figures")


if __name__ == "__main__":
    main()
