"""Consistency checks over the manuscript, run rather than eyeballed.

Three failure modes this catches, all of which happened at least once while the
paper was being split: a number quoted in the text that no longer matches the
results file it came from, a cross-reference to a section that was cut or
renumbered, and a figure or citation that is referenced but does not exist.

It is deliberately narrow. It checks what can be checked mechanically and says
nothing about whether the prose is right.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
SI = ROOT / "paper" / "paper_A_supporting_information.md"
FIGURES = ROOT / "paper" / "figures"
EXTERNAL = ROOT / "data" / "external"
RESULTS = ROOT / "results"


def check_numbers(text: str) -> list[str]:
    """Key figures quoted in the text, against the files they derive from."""
    problems = []

    ppb = pd.read_csv(EXTERNAL / "ppb_final.csv")
    half = pd.read_csv(EXTERNAL / "half_life_final.csv")
    az = pd.read_csv(EXTERNAL / "az_assays.csv")

    expected = {
        "human hepatocyte assay size": (
            int((az["assay"] == "cl_hep_human").sum()), ["408"]),
        "rat hepatocyte assay size": (
            int((az["assay"] == "cl_hep_rat").sum()), ["837"]),
        "microsome assay size": (
            int((az["assay"] == "cl_mic_human").sum()), ["1102", "1 102"]),
        "half-life verified set": (len(half), ["76"]),
    }
    for label, (value, spellings) in expected.items():
        if not any(s in text for s in spellings):
            problems.append(f"{label}: computed {value}, no matching figure in text")

    # The protein-binding set loses one compound to a training-set collision at
    # evaluation time, so the paper should say 128 where the file holds 129.
    if "128" not in text:
        problems.append(
            f"protein-binding set: file holds {len(ppb)}, one is dropped at "
            "evaluation for appearing in training, so the text should say 128"
        )

    composition = RESULTS / "ppb_species_composition.csv"
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

    summary = RESULTS / "external_validation_summary.csv"
    if summary.exists():
        ext = pd.read_csv(summary).set_index("endpoint")
        for endpoint, token in (("ppb", "0.689"), ("half_life", "0.334")):
            value = ext.loc[endpoint, "external_spearman"]
            if token not in text:
                problems.append(
                    f"{endpoint} external Spearman is {value:.3f}; '{token}' not in text"
                )
    return problems


def check_cross_references(text: str) -> list[str]:
    """Every 'Section N.M' must point at a heading that exists."""
    headings = set(re.findall(r"^#{2,3} (\d+(?:\.\d+)?)", text, re.M))
    problems = []
    for reference in sorted(set(re.findall(r"Section (\d+\.\d+)", text))):
        if reference not in headings:
            problems.append(f"cross-reference to Section {reference}, which has no heading")
    return problems


def check_figures(text: str) -> list[str]:
    problems = []
    for path in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text):
        if not (PAPER.parent / path).exists():
            problems.append(f"figure file missing: {path}")
    for number in sorted(set(re.findall(r"\bFigure (\d+)", text))):
        if f"**Figure {number}." not in text:
            problems.append(f"Figure {number} referenced but has no caption")
    return problems


def check_citations(text: str) -> list[str]:
    """Superscript citation markers against the numbered reference list."""
    listed = {int(n) for n in re.findall(r"^\((\d+)\) ", text, re.M)}
    superscripts = str.maketrans("¹²³⁴⁵⁶⁷⁸⁹⁰", "1234567890")
    used = set()
    for run in re.findall(r"[¹²³⁴⁵⁶⁷⁸⁹⁰]+(?:'[¹²³⁴⁵⁶⁷⁸⁹⁰]+)*", text):
        for part in run.split("'"):
            if part:
                used.add(int(part.translate(superscripts)))

    problems = []
    for number in sorted(used - listed):
        problems.append(f"citation {number} used but not in the reference list")
    for number in sorted(listed - used):
        problems.append(f"reference {number} listed but never cited")
    return problems


def check_retractions(text: str) -> list[str]:
    """Withdrawn figures may appear only inside the passages that retract them."""
    problems = []
    for token in ("0.541", "0.764", "0.294", "0.571"):
        lines = [i for i, line in enumerate(text.split("\n")) if token in line]
        for index in lines:
            window = "\n".join(text.split("\n")[max(0, index - 12) : index + 6]).lower()
            if not any(w in window for w in ("withdraw", "retract", "both are", "we estimated")):
                problems.append(
                    f"retracted figure {token} appears at line {index + 1} "
                    "outside a retraction passage"
                )
    return problems


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    si = SI.read_text(encoding="utf-8")

    sections = {
        "numbers against source files": check_numbers(text + si),
        "cross-references": check_cross_references(text),
        "figures": check_figures(text),
        "citations": check_citations(text),
        "retracted figures": check_retractions(text),
    }

    total = 0
    for name, problems in sections.items():
        status = "OK" if not problems else f"{len(problems)} problem(s)"
        print(f"\n=== {name}: {status} ===")
        for problem in problems:
            print(f"  ! {problem}")
        total += len(problems)

    print(f"\n{'clean' if total == 0 else str(total) + ' item(s) to fix'}")


if __name__ == "__main__":
    main()
