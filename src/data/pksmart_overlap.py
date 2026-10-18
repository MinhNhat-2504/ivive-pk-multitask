"""Deduplicate the temporal external set against PKSmart's released test sets.

PKSmart (Seal et al., J Cheminform 2025) publishes two external sets that draw
on FDA approvals, so a compound of ours that also sits in theirs is not
external to a model trained or selected on their data. Checking it is cheap and
the claim does not survive review without it.

Two matching routes, because the two files differ. `FDA_v2.csv` carries drug
names and neutral SMILES; `External_test_315.csv` carries only `smiles_r`, which
is standardised and explicitly protonated -- chlordiazepoxide appears there as
`C/[NH+]=C1/C[NH+](O)...` -- so structures must be neutralised before keying or
every match is missed.

The check also quantifies the apparent-parameter contamination in those files,
since our Section 6 cites it: a true intravenous steady-state volume is almost
never above ~10 L/kg, and a systemic clearance cannot exceed hepatic blood flow
for a compound that is not renally cleared.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"
PKSMART_DIR = Path(
    r"C:\Users\MINHNH~1\AppData\Local\Temp\claude"
    r"\d--Project-AI-Engineer-ivive-pk-multitask"
    r"\2b9907bc-06c9-44fb-aff6-5155ca40a9c3\scratchpad\dl2"
)

#: Above this a reported steady-state volume is almost certainly Vz/F divided
#: by a body weight rather than a true intravenous VDss.
IMPLAUSIBLE_VDSS = 10.0
#: Human hepatic blood flow, mL/min/kg. A systemic clearance above it is not
#: achievable by hepatic extraction alone.
HEPATIC_BLOOD_FLOW = 20.7

_uncharger = rdMolStandardize.Uncharger()


def standard_key(smiles: str | None) -> str | None:
    """Neutralised, salt-stripped InChIKey, so protonation cannot hide a match."""
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if fragments:
        mol = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
    try:
        mol = _uncharger.uncharge(mol)
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    ours = pd.read_csv(EXTERNAL_ROOT / "external_candidates.csv")
    ours = ours[~ours["in_any_tdc_set"]].copy()
    ours["key"] = ours["parent_smiles"].map(standard_key)
    ours = ours.dropna(subset=["key"])

    fda = pd.read_csv(PKSMART_DIR / "FDA_v2.csv")
    bulk = pd.read_csv(PKSMART_DIR / "External_test_315.csv")
    fda["key"] = fda["smiles"].map(standard_key)
    bulk["key"] = bulk["smiles_r"].map(standard_key)

    fda_keys = set(fda["key"].dropna())
    bulk_keys = set(bulk["key"].dropna())
    names = {str(n).strip().lower() for n in fda["name"].dropna()}

    ours["in_pksmart_fda"] = ours["key"].isin(fda_keys)
    ours["in_pksmart_315"] = ours["key"].isin(bulk_keys)
    ours["name_matches_pksmart"] = ours["primary_ingredient"].str.lower().isin(names)
    ours["in_pksmart"] = (
        ours["in_pksmart_fda"] | ours["in_pksmart_315"] | ours["name_matches_pksmart"]
    )

    print(f"\nour temporally clean candidates : {len(ours)}")
    print(f"PKSmart FDA_v2 rows             : {len(fda)} ({fda['key'].notna().sum()} keyed)")
    print(f"PKSmart External_test_315 rows  : {len(bulk)} ({bulk['key'].notna().sum()} keyed)")
    print("\n=== Overlap ===")
    print(f"  structure in FDA_v2           : {int(ours['in_pksmart_fda'].sum())}")
    print(f"  structure in External_test_315: {int(ours['in_pksmart_315'].sum())}")
    print(f"  name matches an FDA_v2 drug   : {int(ours['name_matches_pksmart'].sum())}")
    print(f"  ANY of the three              : {int(ours['in_pksmart'].sum())}")

    shared = ours.loc[ours["in_pksmart"], "primary_ingredient"].tolist()
    if shared:
        print("\n  shared compounds:", ", ".join(sorted(shared)))

    print("\n=== Apparent-parameter contamination in PKSmart's FDA rows ===")
    vdss = pd.to_numeric(fda["Volume_of_distribution_(VDss)"], errors="coerce")
    clear = pd.to_numeric(fda["Clearance_(CL)"], errors="coerce")
    bad_v = fda.loc[vdss > IMPLAUSIBLE_VDSS, ["name", "source"]].assign(vdss=vdss[vdss > IMPLAUSIBLE_VDSS])
    bad_c = fda.loc[clear > HEPATIC_BLOOD_FLOW, ["name", "source"]].assign(cl=clear[clear > HEPATIC_BLOOD_FLOW])
    print(f"  VDss > {IMPLAUSIBLE_VDSS:.0f} L/kg      : {len(bad_v)} of {int(vdss.notna().sum())} "
          f"({100 * len(bad_v) / max(vdss.notna().sum(), 1):.1f}%)")
    print(f"  CL > hepatic blood flow: {len(bad_c)} of {int(clear.notna().sum())} "
          f"({100 * len(bad_c) / max(clear.notna().sum(), 1):.1f}%)")
    print(f"  all implausible rows are source=FDA: "
          f"{set(bad_v['source']) | set(bad_c['source']) == {'FDA'}}")
    print("\n" + bad_v.sort_values("vdss", ascending=False).to_string(index=False))

    out = EXTERNAL_ROOT / "external_candidates_dedup.csv"
    ours.to_csv(out, index=False)
    print(f"\nVERIFIED EXTERNAL, PKSmart-disjoint: {int((~ours['in_pksmart']).sum())}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
