Follow-up comment for issue #405 on `mims-harvard/TDC`.
Post as a new comment on the existing issue, not a new issue.

---

Following up on this with a third observation, on `PPBR_AZ`, found while
checking the first two. It is the same class of issue as the hepatocyte one.

### `PPBR_AZ` is the union of five species' assays

The source deposition, ChEMBL document `CHEMBL3301361`, contains five plasma
protein binding assays rather than one:

| Assay | Species | Records |
|---|---|---|
| `CHEMBL3301365` | human | 1614 |
| `CHEMBL3301366` | rat | 717 |
| `CHEMBL3301367` | dog | 244 |
| `CHEMBL3301368` | mouse | 162 |
| `CHEMBL3301369` | guinea pig | 91 |
| | **total** | **2828 records over 1797 molecules** |

`admet_group` ships `PPBR_AZ` as 2790 rows over exactly 1797 distinct molecules.
Matching each benchmark row against all five assays by (canonical SMILES, value),
97.5% match exactly one assay, and those divide:

| Species | Rows | Share |
|---|---|---|
| human | 1567 | 57.6% |
| rat | 680 | 25.0% |
| dog | 233 | 8.6% |
| mouse | 153 | 5.6% |
| guinea pig | 86 | 3.2% |

So **42.4% of the endpoint is not human plasma**. The dataset page describes it
only as "plasma protein binding rate (AstraZeneca)" and does not mention species.

### Two consequences worth separating

**The documented size does not match the shipped size.** `tdc/metadata.py` has
`"ppbr_az": 1797` and the benchmark page states 1,797 compounds, but
`admet_group.get("PPBR_AZ")` returns 2231 + 559 = 2790 rows. `ADME(name="PPBR_AZ")`
returns 1614. The three numbers are the molecule count, the row count including
the per-species duplicates, and the human assay alone.

**The multiplicity puts a floor under the leaderboard metric.** Because the same
molecule appears up to five times with different values, the test partition is
559 rows over 343 molecules, and 140 of those molecules carry labels that
disagree — by up to 49.3 percentage points. No single prediction can satisfy two
contradictory labels, so taking each molecule's median label as the best possible
prediction gives an irreducible MAE of **1.964 percentage points**. The top ten
entries on the leaderboard are separated from one another by 1.24 points in
total, which is less than that floor.

To be clear about what this is not: **no molecule appears in both the train and
test partitions**, so this is not a leakage report. The split is molecule-
disjoint and the duplication is within partitions.

### Reproducing it

```python
from tdc.benchmark_group import admet_group
g = admet_group(path="data/")
b = g.get("PPBR_AZ")
rows = len(b["train_val"]) + len(b["test"])          # 2790
mols = len(set(b["train_val"].Drug) | set(b["test"].Drug))   # 1797
```

and the five assays are one API call each:

```
https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_chembl_id=CHEMBL3301366
```

Full matching script: https://github.com/MinhNhat-2504/ivive-pk-multitask
(`src/data/ppb_species.py`).

### What would resolve it

Same as for the hepatocyte endpoint — a `species` column, or a documented
human-only variant. For this endpoint there is a second, cheaper option: state
in the dataset description that rows are per measurement rather than per
molecule, so that users and the leaderboard metric are not silently averaging
over cross-species disagreement.

We are happy to open a PR adding species to both endpoints if that would be
useful. As before, if any of this is intentional we will correct our description
of it accordingly.
