"""Renumber paper A's figure references and attach captions.

The figure placeholders were inherited from the combined draft and had gone
stale: they cited Figures 1, 2 and 4, with no Figure 3, because the figure that
had been Figure 3 belongs to the method paper and the one Section 4.5 needs did
not exist yet. This rewrites them as a contiguous set of three with captions.
"""

from __future__ import annotations

from pathlib import Path

PAPER = Path(__file__).resolve().parents[2] / "paper" / "paper_A_benchmark_defects.md"

CAPTIONS = {
    "`[Figure 1]`": (
        "![Figure 1](figures/paperA_fig1_species.png)\n\n"
        "**Figure 1. `Clearance_Hepatocyte_AZ` is built from two species.** Left: "
        "the source deposition holds 408 human and 837 rat hepatocyte clearance "
        "records, each labelled by species in ChEMBL and PubChem. Right: the "
        "1 213 records of the benchmark endpoint, matched against both assays by "
        "value. Two thirds carry the rat value. The benchmark's documentation "
        "separates its clearance data by experiment type and does not mention "
        "species."
    ),
    "`[Figure 2]`": (
        "![Figure 2](figures/paperA_fig2_censoring.png)\n\n"
        "**Figure 2. A quarter to a third of the clearance labels sit on a "
        "detection limit.** Distribution of intrinsic clearance for the two AZ "
        "assays, log scale, with the records at the assay's declared bounds of 3 "
        "and 150 µL/min drawn in colour and the interior in grey. The deposition "
        "states the range as \"<3 to >150\" in its own assay description; ChEMBL "
        "records which values are at a limit in `standard_relation` and the "
        "benchmark's conversion drops that column. The next most frequent label "
        "in either endpoint occurs three to six times."
    ),
    "`[Figure 4]`": (
        "![Figure 3](figures/paperA_fig3_calibration.png)\n\n"
        "**Figure 3. Training on censored labels compresses predictions, and the "
        "cost falls on the unbound fraction.** Left: predicted against observed "
        "plasma protein binding for the 128 external compounds, with the assay's "
        "declared bounds as vertical lines and the compounds reported at a bound "
        "as open circles. Weakly bound compounds are over-predicted by up to 36 "
        "percentage points and highly bound ones under-predicted by about 6; the "
        "model lands on the correct side of a bound in 4% of right-censored "
        "cases. Right: the same predictions as unbound fraction, "
        "$f_u = 1 - \\mathrm{PPB}/100$, for the 37 compounds bound above 99%. "
        "Because $f_u$ is a difference from 100, a few percentage points of "
        "compression is a median error of 8.3-fold in the quantity that every "
        "well-stirred clearance calculation consumes. Dashed line: identity. "
        "Dotted: ten-fold."
    ),
}


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    for placeholder, caption in CAPTIONS.items():
        if placeholder not in text:
            raise SystemExit(f"placeholder missing: {placeholder}")
        text = text.replace(placeholder, caption)
    PAPER.write_text(text, encoding="utf-8")
    print(f"attached {len(CAPTIONS)} captions; figures renumbered 1-3")


if __name__ == "__main__":
    main()
