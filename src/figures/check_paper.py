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
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"
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


def check_numbers_paper_c(text: str) -> list[str]:
    """Paper C's quoted figures, against the results files they come from."""
    problems = []

    def require(label: str, value: float, spelling: str) -> None:
        if spelling not in text:
            problems.append(f"{label}: computed {value}, '{spelling}' not in text")

    compression = RESULTS / "ppb_compression_summary.csv"
    if compression.exists():
        table = pd.read_csv(compression).set_index("model")
        for model, fold, over in (
            ("ridge (linear)", "20.0", "8.2"),
            ("SVR (RBF kernel)", "11.2", "11.2"),
            ("k-NN (k=5)", "17.9", "17.9"),
            ("random forest", "10.1", "10.1"),
            ("gradient boosting", "8.3", "8.3"),
            ("neural network", "26.5", "26.5"),
        ):
            require(f"{model} fu fold error",
                    round(float(table.loc[model, "fu_fold_error_median"]), 1), fold)
            require(f"{model} fu over-estimate",
                    round(float(table.loc[model, "fu_overprediction_median"]), 1), over)
        for model, spelling in (("ridge (linear)", "| **14** |"),
                                ("neural network", "| **6**  |")):
            require(f"{model} predictions at the assay ceiling",
                    int(table.loc[model, "n_at_assay_ceiling"]), spelling)
        for model, spelling in (("gradient boosting", "35 / 37"),
                                ("ridge (linear)", "23 / 37"),
                                ("neural network", "31 / 37")):
            require(f"{model} over-predicted count",
                    int(table.loc[model, "n_overpredicted"]), spelling)

    scale = RESULTS / "target_scale_per_seed.csv"
    if scale.exists():
        frame = pd.read_csv(scale).groupby("scale")
        for column, which, spelling in (
            ("fu_fold_median_high", "percent", "8.71"),
            ("fu_fold_median_high", "log_fu", "2.35"),
            ("fu_fold_median_high", "logit", "2.46"),
            ("frac_within_2fold_high", "percent", "0.077"),
            ("frac_within_2fold_high", "log_fu", "0.408"),
            ("mae_percent", "logit", "10.15"),
        ):
            require(f"target scale {which} {column}",
                    round(float(frame.get_group(which)[column].mean()), 3), spelling)

    families = RESULTS / "scale_across_families_per_seed.csv"
    if families.exists():
        table = pd.read_csv(families).pivot_table(
            index="family", columns="scale", values="fu_fold_median_high"
        )
        for family, spelling in (
            ("ridge (linear)", "20.00"), ("SVR (RBF kernel)", "15.32"),
            ("k-NN (k=5)", "17.90"), ("random forest", "10.50"),
            ("gradient boosting", "8.28"), ("neural network", "17.38"),
        ):
            require(f"{family} percent-scale fold error",
                    round(float(table.loc[family, "percent"]), 2), spelling)

    submission = RESULTS / "leaderboard_submission_per_seed.csv"
    if submission.exists():
        frame = pd.read_csv(submission)
        official = frame[frame["seed"].isin([1, 2, 3, 4, 5])].groupby("scale")["value"]
        for which, spelling in (("percent", "8.009"), ("log_fu", "7.698"),
                                ("logit", "7.553")):
            require(f"official-protocol MAE, {which}",
                    round(float(official.get_group(which).mean()), 3), spelling)

    human = RESULTS / "human_only_scale_per_seed.csv"
    if human.exists():
        frame = pd.read_csv(human)
        frame = frame[frame["view"] == "external"].groupby(["training", "scale"])
        for training, scale_name, spelling in (
            ("pooled (as distributed)", "percent", "8.710"),
            ("human only", "percent", "8.779"),
            ("human only", "log_fu", "2.158"),
            ("human only", "logit", "2.154"),
        ):
            value = frame.get_group((training, scale_name))["fu_fold_median_high"].mean()
            require(f"{training} / {scale_name}", round(float(value), 3), spelling)

    tuned = RESULTS / "tuned_scale_per_seed.csv"
    if tuned.exists():
        frame = pd.read_csv(tuned)
        for scale_name, state, column, spelling in (
            ("percent", "tuned", "test_mae", "7.991"),
            ("logit", "tuned", "test_mae", "7.410"),
            ("log_fu", "tuned", "test_mae", "7.611"),
            ("percent", "tuned", "ext_fu_fold_median_high", "8.70"),
            ("percent", "untuned", "ext_fu_fold_median_high", "8.97"),
            ("logit", "untuned", "ext_fu_fold_median_high", "2.59"),
            ("log_fu", "untuned", "ext_fu_fold_median_high", "2.45"),
            ("log_fu", "tuned", "ext_fu_fold_median_high", "2.65"),
            ("log_fu", "tuned", "ext_frac_within_2fold_high", "0.330"),
            ("log_fu", "untuned", "ext_frac_within_2fold_high", "0.384"),
        ):
            block = frame[(frame["scale"] == scale_name)
                          & (frame["tuning"] == state)][column]
            require(f"{scale_name} {state} {column}",
                    round(float(block.mean()), 4), spelling)

    propagation = RESULTS / "fu_propagation_parametric.csv"
    if propagation.exists():
        frame = pd.read_csv(propagation)
        frame = frame[frame["fu_fold_error"] == 8.3]
        for fu, column, spelling in (
            (0.001, "CLh_fold_error_at_CLint_150", "7.96"),
            (0.050, "CLh_fold_error_at_CLint_150", "3.10"),
            (0.001, "CLh_fold_error_at_CLint_3", "8.29"),
        ):
            value = float(frame.loc[frame["fu_true"] == fu, column].iloc[0])
            require(f"propagation at fu={fu}, {column}", round(value, 2), spelling)

    composition = RESULTS / "ppb_species_composition.csv"
    if composition.exists():
        species = pd.read_csv(composition).set_index("species")
        for name, spelling in (("human", "1567"), ("rat", "680"), ("dog", "233"),
                               ("mouse", "153"), ("guinea pig", "86")):
            require(f"{name} unambiguous rows",
                    int(species.loc[name, "unambiguous_rows"]), spelling)

    floor_file = RESULTS / "ppb_metric_floor.csv"
    if floor_file.exists():
        floor = pd.read_csv(floor_file).iloc[0]
        require("MAE floor", floor["mae_floor"], "1.964")
        require("RMSE floor", floor["rmse_floor"], "4.478")

    return problems


def check_cross_references(text: str) -> list[str]:
    """Every 'Section N.M' must point at a heading that exists.

    A section of an external document -- the FDA label's pharmacokinetics
    section, say -- is written with the section sign, as §12.3, so that this
    check can stay strict about the spelled-out form.
    """
    headings = set(re.findall(r"^#{2,3} (\d+(?:\.\d+)?)", text, re.M))
    problems = []
    for reference in sorted(set(re.findall(r"Section (\d+\.\d+)", text))):
        if reference not in headings:
            problems.append(f"cross-reference to Section {reference}, which has no heading")
    return problems


def check_figures(text: str) -> list[str]:
    problems = []
    for path in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text):
        if not (ROOT / "paper" / path).exists():
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
    paper = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PAPER
    if not paper.is_absolute():
        paper = ROOT / paper
    text = paper.read_text(encoding="utf-8")
    is_paper_a = paper.name == DEFAULT_PAPER.name

    sections = {
        "cross-references": check_cross_references(text),
        "figures": check_figures(text),
        "citations": check_citations(text),
    }
    if not is_paper_a:
        # Paper C's supporting information carries the full tables, so a figure
        # may be quoted there and nowhere else; the checks run over both.
        si_c = ROOT / "paper" / "paper_C_supporting_information.md"
        combined = text + (si_c.read_text(encoding="utf-8") if si_c.exists() else "")
        sections = {
            "numbers against source files": check_numbers_paper_c(combined),
            **sections,
        }
    if is_paper_a:
        # The quoted figures and the retraction passages are paper A's, and the
        # supporting information belongs to it.
        si = SI.read_text(encoding="utf-8")
        sections = {
            "numbers against source files": check_numbers(text + si),
            **sections,
            "retracted figures": check_retractions(text),
        }
    print(f"checking {paper.name}")

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
