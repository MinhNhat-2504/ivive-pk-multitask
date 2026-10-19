"""Correct paper A's claim about PPBR_AZ, and add the second species pooling.

Section 4.1 asserted that `Clearance_Microsome_AZ` and `PPBR_AZ` "match a single
human assay exactly, at 1 102 and 1 614 records". That is right for the
microsome endpoint and wrong for protein binding: 1 614 is the size of the
upstream human assay, not of the matched benchmark subset, and `PPBR_AZ` holds
2 790 rows over 1 797 molecules. The deposition contains five protein binding
assays, one per species, and the benchmark is their union.

That makes species pooling a finding about two endpoints rather than one, so it
extends section 4.1 rather than becoming a third defect; the paper's title and
its count of two defects stand. The within-molecule duplication the pooling
creates also puts a floor under the leaderboard metric, which belongs in 4.3
where the precision of the evaluation is discussed.

Numbers are from `results/ppb_all_species_upstream.csv` and the value match in
`notes/22`. `check_paper.py` verifies them against those files.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"

CORRECTION = (
    """Species is not mentioned. `Clearance_Microsome_AZ` and
`PPBR_AZ` match a single human assay exactly, at 1 102 and 1 614 records.""",
    """Species is not mentioned. `Clearance_Microsome_AZ` does match a single
human assay exactly, at 1 102 records.

**The same pooling occurs in `PPBR_AZ`, across five species.** The deposition
contains five plasma protein binding assays, one per species: human, rat, dog,
mouse and guinea pig, holding 1 614, 717, 244, 162 and 91 records over 1 797
distinct molecules. `PPBR_AZ` as distributed by the benchmark group holds 2 790
rows over exactly 1 797 distinct molecules, and the benchmark's own metadata
records the dataset size as 1 797. Matching each row against all five assays by
value, 97.5% match exactly one, and those divide as

| | records | share |
|---|---|---|
| Human plasma | 1 567 | 57.6% |
| **Rat plasma** | **680** | **25.0%** |
| **Dog plasma** | **233** | **8.6%** |
| **Mouse plasma** | **153** | **5.6%** |
| **Guinea pig plasma** | **86** | **3.2%** |

so 42.4% of the endpoint is non-human. The upstream human assay holds one record
per molecule with no duplication at all; the 691 molecules that appear two to
five times in the benchmark do so because the same compound was measured in
several species. No molecule crosses the train/test boundary, so this is not a
leakage finding.

The two endpoints differ in one respect worth stating. For hepatocyte clearance
the pooled assays are an in vitro measurement in two species; for protein
binding they are five. In both cases the species label exists upstream, is
recorded in ChEMBL, and is absent from the distributed benchmark.""",
)

FLOOR = (
    """### 4.3 A reconciliation, and the precision the evaluation supports""",
    """### 4.3 A reconciliation, and the precision the evaluation supports

The species pooling in `PPBR_AZ` has a direct consequence for its leaderboard.
Its test partition is 559 rows over 343 molecules, and 140 of those molecules
carry two or more labels that disagree, by up to 49.3 percentage points. No
single prediction can satisfy two contradictory labels, so the row-weighted mean
absolute error has a floor: taking each molecule's median label as the best
possible prediction, **1.964 percentage points of the metric cannot be removed by
any model** (4.478 for root mean squared error, using the mean). Published
entries on this endpoint sit near 7.4 to 7.8 percentage points and are separated
from one another by a few tenths, so the irreducible component is several times
the spread being competed over, and it is cross-species disagreement rather than
assay noise.""",
)


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    for old, new in (CORRECTION, FLOOR):
        if old not in text:
            raise SystemExit(f"anchor not found:\n  {old[:90]}")
        if text.count(old) != 1:
            raise SystemExit(f"anchor is not unique:\n  {old[:90]}")
        text = text.replace(old, new, 1)
    PAPER.write_text(text, encoding="utf-8")
    print("paper A: PPBR_AZ claim corrected, five-species pooling and MAE floor added")


if __name__ == "__main__":
    main()
