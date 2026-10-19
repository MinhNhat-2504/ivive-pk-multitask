"""Register paper C's quoted figures with the consistency checker.

Paper C was getting the structural checks -- cross-references, figures,
citations -- and none of the number checks, because those were written against
paper A's results files. Every figure corrected in the last two days was the
kind this catches: a statistic whose definition changed in the code while the
manuscript kept the old value.

Each entry names the results file a figure comes from, recomputes it, and
requires the spelling to appear in the text. Presence is a weak test, but it is
the one that fails when a number drifts: recompute, and the new value is not in
the manuscript.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "src" / "figures" / "check_paper.py"

FUNCTION = '''

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
            ("ridge (linear)", "30.0", "8.2"),
            ("SVR (RBF kernel)", "11.2", "11.2"),
            ("k-NN (k=5)", "17.9", "17.9"),
            ("random forest", "10.1", "10.1"),
            ("gradient boosting", "8.8", "8.3"),
            ("neural network", "36.3", "26.5"),
        ):
            require(f"{model} fu fold error",
                    round(float(table.loc[model, "fu_fold_error_median"]), 1), fold)
            require(f"{model} fu over-estimate",
                    round(float(table.loc[model, "fu_overprediction_median"]), 1), over)
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
'''

ANCHOR = """def check_cross_references(text: str) -> list[str]:"""

WIRING = ('''    if is_paper_a:''',
          '''    if not is_paper_a:
        sections = {
            "numbers against source files": check_numbers_paper_c(text),
            **sections,
        }
    if is_paper_a:''')


def main() -> None:
    text = CHECKER.read_text(encoding="utf-8")
    if "check_numbers_paper_c" in text:
        print("already registered")
        return
    if text.count(ANCHOR) != 1 or text.count(WIRING[0]) != 1:
        raise SystemExit("anchors not found exactly once in check_paper.py")
    text = text.replace(ANCHOR, FUNCTION.strip() + "\n\n\n" + ANCHOR, 1)
    text = text.replace(WIRING[0], WIRING[1], 1)
    CHECKER.write_text(text, encoding="utf-8")
    print("paper C's figures registered with the checker")


if __name__ == "__main__":
    main()
