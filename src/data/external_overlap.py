"""Does the temporal external set actually sit outside the training data?

The whole value of a post-2017 approval set is that no model benchmarked on
these leaderboards can have seen it. That has to be demonstrated rather than
argued from approval dates: a compound approved in 2019 may have been in the
literature, and in ChEMBL, for a decade beforehand.

Three checks, increasingly strict:

    InChIKey          the same compound, by connectivity and stereochemistry
    InChIKey skeleton the same compound ignoring stereochemistry and salt form
    scaffold          the same Bemis-Murcko framework

Drug-likeness filtering happens here too, since the FDA set contains peptides,
oligonucleotides, radiopharmaceuticals and imaging agents that no small
molecule model should be asked to predict.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from src.data.load_tdc import PK_ENDPOINTS, load_task

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"

#: Drug-likeness bounds. Deliberately loose -- the intent is to remove
#: peptides, oligonucleotides and contrast agents, not to impose Lipinski on
#: compounds that were approved despite it.
MAX_MW = 900.0
MIN_MW = 100.0
ALLOWED_ELEMENTS = set("C H N O S P F Cl Br I Na K Ca Mg".split())


def largest_fragment(smiles: str) -> str | None:
    """Strip counterions by keeping the heaviest fragment."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not fragments:
        return None
    heaviest = max(fragments, key=Descriptors.MolWt)
    return Chem.MolToSmiles(heaviest)


def is_small_molecule(smiles: str) -> tuple[bool, str]:
    """Drug-like organic small molecule, with the reason when it is not."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "unparsable"
    mw = Descriptors.MolWt(mol)
    if mw > MAX_MW:
        return False, f"mw {mw:.0f} > {MAX_MW:.0f}"
    if mw < MIN_MW:
        return False, f"mw {mw:.0f} < {MIN_MW:.0f}"
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    exotic = elements - ALLOWED_ELEMENTS
    if exotic:
        return False, f"elements {sorted(exotic)}"
    # Peptides slip past the mass bound when short; an amide-rich backbone with
    # many rotatable bonds and few rings is the usual signature.
    if Chem.rdMolDescriptors.CalcNumAmideBonds(mol) >= 5 and mol.GetRingInfo().NumRings() <= 2:
        return False, "peptide-like"
    return True, ""


def keys_for(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    inchikey = Chem.MolToInchiKey(mol)
    try:
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        scaffold = ""
    return {
        "inchikey": inchikey,
        "inchikey_skeleton": inchikey.split("-")[0] if inchikey else "",
        "scaffold": scaffold,
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    frame = pd.read_csv(EXTERNAL_ROOT / "nme_labels_structures.csv")
    frame = frame[frame["smiles"].notna() & (frame["smiles"] != "")].copy()
    log.info("%d applications with a structure", len(frame))

    frame["parent_smiles"] = frame["smiles"].map(largest_fragment)
    frame = frame.dropna(subset=["parent_smiles"])

    verdicts = frame["parent_smiles"].map(is_small_molecule)
    frame["is_small_molecule"] = [v[0] for v in verdicts]
    frame["exclusion_reason"] = [v[1] for v in verdicts]

    print("\n=== Drug-likeness filter ===")
    print(f"with a structure      : {len(frame)}")
    print(f"small molecules       : {int(frame['is_small_molecule'].sum())}")
    rejected = frame.loc[~frame["is_small_molecule"], "exclusion_reason"]
    print("excluded, by reason:")
    print(rejected.str.replace(r"mw .*", "mass out of range", regex=True).value_counts().to_string())

    small = frame[frame["is_small_molecule"]].copy()
    keys = small["parent_smiles"].map(keys_for).apply(pd.Series)
    small = pd.concat([small.reset_index(drop=True), keys.reset_index(drop=True)], axis=1)
    small = small[small["inchikey"].notna() & (small["inchikey"] != "")]
    small = small.drop_duplicates(subset="inchikey")
    print(f"\nunique small molecules: {len(small)}")

    print("\n=== Overlap with TDC training data ===")
    rows = []
    for tdc_name, key in PK_ENDPOINTS.items():
        task = load_task(tdc_name, key)
        tdc_keys = task.frame["smiles"].map(keys_for).apply(pd.Series)
        row = {"endpoint": key, "n_tdc": len(task)}
        for level in ("inchikey", "inchikey_skeleton", "scaffold"):
            reference = set(tdc_keys[level]) - {""}
            overlap = small[level].isin(reference)
            row[f"shared_{level}"] = int(overlap.sum())
        rows.append(row)
    overlap_table = pd.DataFrame(rows)
    print(overlap_table.to_string(index=False))

    any_exact = pd.Series(False, index=small.index)
    for tdc_name, key in PK_ENDPOINTS.items():
        task = load_task(tdc_name, key)
        tdc_keys = task.frame["smiles"].map(keys_for).apply(pd.Series)
        any_exact |= small["inchikey_skeleton"].isin(set(tdc_keys["inchikey_skeleton"]) - {""})

    small["in_any_tdc_set"] = any_exact
    clean = small[~small["in_any_tdc_set"]]
    print(f"\nin at least one TDC set (skeleton match): {int(any_exact.sum())}")
    print(f"TEMPORALLY CLEAN external candidates     : {len(clean)}")
    print(f"  of those with a PK section             : {int(clean['has_pk_section'].sum())}")

    out = EXTERNAL_ROOT / "external_candidates.csv"
    small.to_csv(out, index=False)
    overlap_table.to_csv(EXTERNAL_ROOT / "external_overlap.csv", index=False)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
