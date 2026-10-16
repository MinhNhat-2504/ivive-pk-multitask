"""Load the five human PK endpoints of TDC ADMET and report task overlap.

The five endpoints are the ones coupled by the well-stirred IVIVE relation:

    CL_hep, CL_mic   in vitro intrinsic clearance  (uL/min/10^6 cells, uL/min/mg)
    PPBR             fraction bound to plasma protein  -> f_u = 1 - PPBR/100
    VDss             steady-state volume of distribution (L/kg)
    Half_Life        terminal half-life (h)        t_1/2 = ln(2) * VDss / CL

Overlap matters: the IVIVE constraint can only be applied to molecules that
carry labels for several of these at once. This module reports that overlap so
the constraint's coverage is a measured number, not an assumption.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from tdc.single_pred import ADME

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

#: TDC dataset name -> short task key used throughout the project.
PK_ENDPOINTS: dict[str, str] = {
    "Clearance_Hepatocyte_AZ": "cl_hep",
    "Clearance_Microsome_AZ": "cl_mic",
    "PPBR_AZ": "ppbr",
    "VDss_Lombardo": "vdss",
    "Half_Life_Obach": "half_life",
}

DATA_ROOT = Path(__file__).resolve().parents[2] / "data"


@dataclass(frozen=True)
class Task:
    """One PK endpoint, keyed by canonical SMILES."""

    key: str
    tdc_name: str
    frame: pd.DataFrame  # columns: smiles, value

    def __len__(self) -> int:
        return len(self.frame)


def canonical_smiles(smiles: str) -> str | None:
    """Return the RDKit canonical SMILES, or None if the string does not parse.

    Canonicalising before the merge is what makes the overlap count real: TDC
    draws these five sets from different primary sources, so the same molecule
    routinely arrives written two different ways.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def load_task(tdc_name: str, key: str, root: Path = DATA_ROOT) -> Task:
    """Download (or reuse the cache of) one TDC endpoint and canonicalise it."""
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    frame = ADME(name=tdc_name, path=str(raw_dir)).get_data()
    frame = frame.rename(columns={"Drug": "smiles", "Y": "value"})[["smiles", "value"]]

    frame["smiles"] = frame["smiles"].map(canonical_smiles)
    unparsed = int(frame["smiles"].isna().sum())
    if unparsed:
        log.warning("%s: dropped %d unparsable SMILES", tdc_name, unparsed)
    frame = frame.dropna(subset=["smiles"])

    # A molecule measured twice in the same assay is averaged rather than
    # dropped: both readings are evidence, and the spread between them is the
    # first handle we have on the assay's own noise floor.
    frame = frame.groupby("smiles", as_index=False)["value"].mean()

    return Task(key=key, tdc_name=tdc_name, frame=frame)


def load_all(root: Path = DATA_ROOT) -> dict[str, Task]:
    """Load every endpoint in :data:`PK_ENDPOINTS`."""
    return {key: load_task(name, key, root) for name, key in PK_ENDPOINTS.items()}


def build_matrix(tasks: dict[str, Task]) -> pd.DataFrame:
    """Outer-join every task into one smiles x endpoint matrix of labels."""
    matrix: pd.DataFrame | None = None
    for key, task in tasks.items():
        column = task.frame.rename(columns={"value": key})
        matrix = column if matrix is None else matrix.merge(column, on="smiles", how="outer")
    assert matrix is not None, "PK_ENDPOINTS must not be empty"
    return matrix


def overlap_report(matrix: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Summarise how much label co-occurrence the IVIVE constraint can use."""
    keys = [c for c in matrix.columns if c != "smiles"]
    present = matrix[keys].notna()

    pairwise = pd.DataFrame(
        [[int((present[a] & present[b]).sum()) for b in keys] for a in keys],
        index=keys,
        columns=keys,
    )

    n_labels = present.sum(axis=1)
    coverage = (
        n_labels.value_counts()
        .sort_index()
        .rename_axis("n_endpoints")
        .reset_index(name="n_molecules")
    )

    # The subsets each physical relation actually needs.
    constraint_rows = pd.DataFrame(
        [
            {
                "constraint": "well-stirred: CL_hep + f_u -> CL_int,in vivo",
                "needs": "cl_hep + ppbr",
                "n_molecules": int((present["cl_hep"] & present["ppbr"]).sum()),
            },
            {
                "constraint": "well-stirred: CL_mic + f_u -> CL_int,in vivo",
                "needs": "cl_mic + ppbr",
                "n_molecules": int((present["cl_mic"] & present["ppbr"]).sum()),
            },
            {
                "constraint": "t_1/2 = ln2 * VDss / CL  (closes the loop)",
                "needs": "half_life + vdss + cl_hep",
                "n_molecules": int(
                    (present["half_life"] & present["vdss"] & present["cl_hep"]).sum()
                ),
            },
            {
                "constraint": "full chain, all five labels",
                "needs": " + ".join(keys),
                "n_molecules": int(present.all(axis=1).sum()),
            },
        ]
    )

    return {"pairwise": pairwise, "coverage": coverage, "constraints": constraint_rows}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    tasks = load_all()
    print("\n=== Per-endpoint sizes (after canonicalisation + dedup) ===")
    for key, task in tasks.items():
        print(f"  {key:10s} {task.tdc_name:26s} n = {len(task):5d}")

    matrix = build_matrix(tasks)
    print(f"\nUnion of all endpoints: {len(matrix)} unique molecules")

    report = overlap_report(matrix)
    print("\n=== Pairwise label co-occurrence ===")
    print(report["pairwise"].to_string())
    print("\n=== Molecules by number of endpoints labelled ===")
    print(report["coverage"].to_string(index=False))
    print("\n=== Coverage of each IVIVE constraint ===")
    print(report["constraints"].to_string(index=False))

    out = DATA_ROOT / "processed"
    out.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out / "pk_matrix.csv", index=False)
    print(f"\nWrote {out / 'pk_matrix.csv'}")


if __name__ == "__main__":
    main()
