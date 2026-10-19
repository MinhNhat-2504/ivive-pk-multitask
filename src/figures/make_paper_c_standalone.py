"""Remove paper C's load-bearing dependency on the companion paper.

Paper C cites the companion paper eight times. Five of those are ordinary
supporting citations -- a methodological justification, a comparison of two
statistics, a remark about censoring. Three are load-bearing: the species
composition of the endpoint, which Sections 3.2 and 3.6 rest on, and the claim
that the microsomal clearance assay is a single species, which the propagation
in Section 3.1 rests on. A reviewer is entitled to object to a claim of that
weight resting on a reference marked "Submitted".

So the derivation is restated in Methods, with the assay identifiers and the
match procedure, and a reader can reproduce it without the companion paper. The
finding is not claimed here: the text says where it is reported in full and why
it is repeated.

The microsome claim is cheaper still -- the assay identifier and the record
count settle it in a clause.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_C_target_scale.md"

SPECIES_BLOCK = """
**Species composition.** The endpoint is not a single assay. Its source
deposition, ChEMBL document `CHEMBL3301361`, contains five plasma protein
binding assays measured by equilibrium dialysis, differing only in the plasma
used: `CHEMBL3301365` human (1614 records), `CHEMBL3301366` rat (717),
`CHEMBL3301367` dog (244), `CHEMBL3301368` mouse (162) and `CHEMBL3301369`
guinea pig (91) — 2828 records over 1797 distinct molecules.[[chembl]]
`PPBR_AZ` holds 2790 rows over exactly 1797 distinct molecules.

Attributing each benchmark row to an assay requires only a match on
(canonical SMILES, value), since the benchmark preserves the deposited values
unchanged. 97.5% of rows match exactly one of the five assays, 1.4% match more
than one because two species happened to give the same value for a molecule,
and 1.2% match none. Of the unambiguous rows:

| Plasma | Rows | Share |
|---|---|---|
| Human | 1567 | 57.6% |
| Rat | 680 | 25.0% |
| Dog | 233 | 8.6% |
| Mouse | 153 | 5.6% |
| Guinea pig | 86 | 3.2% |

so 42.4% of the endpoint is not human. The benchmark's dataset page describes it
as "plasma protein binding rate (AstraZeneca)" and states its size as 1797
compounds; neither species nor the row multiplicity is mentioned.[[leaderboard]]
No molecule appears in both the train and test partitions, so the duplication is
within partitions and is not a leakage finding.

This is reported, with the equivalent analysis for the hepatocyte clearance
endpoint, in a companion paper.[[prior]] It is re-derived here because Sections
3.2 and 3.6 depend on it and the present paper should stand without the other.

### 2.2 Features and model families"""

EDITS = [
    ("""### 2.2 Features and model families""", SPECIES_BLOCK),
    # The intro no longer leans on the companion paper for the composition.
    (
        """benchmark is their union, 42.4% of it non-human.¹⁰ The same molecule""",
        """benchmark is their union, 42.4% of it non-human (Section 2.1). The same molecule""",
    ),
    # Section 3.2 likewise.
    (
        """`PPBR_AZ` is the union of five plasma protein binding assays from its source
deposition — human, rat, dog, mouse and guinea pig — and 42.4% of its rows are
not human.¹⁰ The same molecule consequently appears up to five times with
different values.""",
        """`PPBR_AZ` is the union of five plasma protein binding assays from its source
deposition — human, rat, dog, mouse and guinea pig — and 42.4% of its rows are
not human (Section 2.1). The same molecule consequently appears up to five times
with different values.""",
    ),
    (
        """any model** (4.478 for root mean squared error, using the mean). No molecule
crosses the train/test boundary, so this is not a leakage finding.¹⁰""",
        """any model** (4.478 for root mean squared error, using the mean). No molecule
crosses the train/test boundary, so this is not a leakage finding.""",
    ),
    # The microsome assay's species is settled by naming the assay.
    (
        """$\mathrm{CL_{int}}$ from the human microsomal assay of the same deposition —
which, unlike its hepatocyte counterpart, is a single species¹⁰ — and
combines it in turn with measured and predicted $f_u$.""",
        """$\mathrm{CL_{int}}$ from the human microsomal assay of the same deposition,
`CHEMBL3301370`, whose 1102 records `Clearance_Microsome_AZ` reproduces exactly
and which — unlike the hepatocyte assay of the same deposition[[prior]] — is a
single species, and combines it in turn with measured and predicted $f_u$.""",
    ),
]


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    for old, new in EDITS:
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique or missing:\n  {old[:90]}")
        text = text.replace(old, new, 1)
    PAPER.write_text(text, encoding="utf-8")
    print(f"paper C: {len(EDITS)} passage(s) made self-contained")


if __name__ == "__main__":
    main()
