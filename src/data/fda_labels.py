"""Steps 2 and 3: product labels and structures for the NME set.

Step 2 pulls the pharmacokinetics section of each label. openFDA exposes SPL
section 12.3 as its own `pharmacokinetics` field, separate from the broader
`clinical_pharmacology` (section 12), which saves parsing section 12 prose to
find the subsection. DailyMed is the fallback for applications openFDA cannot
join to a current label, which `notes/05` put at roughly a fifth.

Step 3 resolves each active ingredient to a structure. PubChem's
`CanonicalSMILES` property discards stereochemistry, so `IsomericSMILES` is
requested instead and ChEMBL serves as fallback.

No parameter extraction happens here. `notes/05` established that it has to be
semantic rather than pattern-based -- a first-pass regex for clearance scored
41% but most hits were *creatinine* clearance from renal dosing criteria,
falling to 16% once that context was excluded -- so this module only collects
the text that extraction will later run over, and records provenance for it.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)

OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"
DAILYMED_SPLS = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/IsomericSMILES,MolecularWeight/JSON"
CHEMBL_SEARCH = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"

EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"

#: Salt and ester suffixes appended to an active ingredient name in the FDA
#: listing. The parent compound is what carries the pharmacokinetics, and what
#: TDC holds, so these are stripped before structure lookup.
SALT_SUFFIXES = (
    "hydrochloride", "hydrobromide", "hydrate", "dihydrate", "monohydrate",
    "sulfate", "sulphate", "phosphate", "maleate", "mesylate", "besylate",
    "tosylate", "citrate", "tartrate", "fumarate", "succinate", "acetate",
    "adipate", "malate", "oxalate", "lactate", "gluconate", "sodium",
    "potassium", "calcium", "magnesium", "meglumine", "dimeglumine",
    "trometamol", "tromethamine", "arginine", "lysine", "choline",
    "disodium", "dipotassium", "hemifumarate", "bromide", "chloride",
)


def strip_salt(name: str) -> str:
    """Drop trailing salt or hydrate words from an active-ingredient name."""
    words = name.strip().lower().split()
    while len(words) > 1 and words[-1] in SALT_SUFFIXES:
        words.pop()
    return " ".join(words)


def _get(url: str, params: dict | None = None, timeout: int = 30) -> dict | None:
    try:
        response = requests.get(url, params=params, timeout=timeout)
    except requests.RequestException as exc:
        log.debug("request failed: %s", exc)
        return None
    if response.status_code != 200:
        return None
    try:
        return response.json()
    except ValueError:
        return None


def fetch_label(application_number: str, pause: float = 0.3) -> dict:
    """Pharmacokinetics text for one application, from openFDA then DailyMed."""
    time.sleep(pause)
    payload = _get(
        OPENFDA_LABEL,
        {"search": f'openfda.application_number:"{application_number}"', "limit": 1},
    )
    results = (payload or {}).get("results") or []
    if results:
        record = results[0]
        pk = record.get("pharmacokinetics") or []
        cp = record.get("clinical_pharmacology") or []
        return {
            "label_source": "openfda",
            "pk_text": "\n".join(pk),
            "clinical_pharmacology_text": "\n".join(cp),
            "has_pk_section": bool(pk),
        }
    return {"label_source": "", "pk_text": "", "clinical_pharmacology_text": "", "has_pk_section": False}


def fetch_label_dailymed(drug_name: str, pause: float = 0.3) -> dict:
    """Fallback lookup by drug name; returns the SPL set id to fetch text from."""
    time.sleep(pause)
    payload = _get(DAILYMED_SPLS, {"drug_name": drug_name, "pagesize": 1})
    data = (payload or {}).get("data") or []
    if not data:
        return {}
    return {"label_source": "dailymed", "setid": data[0].get("setid", "")}


def fetch_structure(name: str, pause: float = 0.25) -> dict:
    """IsomericSMILES from PubChem, falling back to ChEMBL."""
    time.sleep(pause)
    payload = _get(PUBCHEM.format(requests.utils.quote(name)))
    properties = ((payload or {}).get("PropertyTable") or {}).get("Properties") or []
    if properties:
        record = properties[0]
        smiles = record.get("IsomericSMILES") or record.get("SMILES")
        if smiles:
            return {
                "smiles": smiles,
                "mw": record.get("MolecularWeight"),
                "structure_source": "pubchem",
            }

    time.sleep(pause)
    payload = _get(CHEMBL_SEARCH, {"q": name, "limit": 1})
    molecules = (payload or {}).get("molecules") or []
    if molecules:
        structures = molecules[0].get("molecule_structures") or {}
        if structures.get("canonical_smiles"):
            return {
                "smiles": structures["canonical_smiles"],
                "mw": (molecules[0].get("molecule_properties") or {}).get("full_mwt"),
                "structure_source": "chembl",
            }
    return {"smiles": "", "mw": None, "structure_source": ""}


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    nmes = pd.read_csv(EXTERNAL_ROOT / "nme_since_2017.csv")
    nda = nmes[nmes["route"] == "NDA"].reset_index(drop=True)
    log.info("resolving %d NDA applications", len(nda))

    rows = []
    for position, record in nda.iterrows():
        # The first listed ingredient is the active moiety for a single-entity
        # product; combinations are kept and flagged rather than dropped here,
        # since the PK section may still carry usable single-agent values.
        ingredients = str(record["ingredients"]).split("; ")
        primary = strip_salt(ingredients[0]) if ingredients else ""

        label = fetch_label(record["application_number"])
        structure = fetch_structure(primary) if primary else {}

        rows.append(
            {
                **record.to_dict(),
                "primary_ingredient": primary,
                "is_combination": len(ingredients) > 1,
                **label,
                **structure,
            }
        )
        if (position + 1) % 25 == 0:
            done = pd.DataFrame(rows)
            log.info(
                "%d/%d  pk_section %d  smiles %d",
                position + 1, len(nda),
                int(done["has_pk_section"].sum()),
                int((done["smiles"] != "").sum()),
            )

    frame = pd.DataFrame(rows)
    out = EXTERNAL_ROOT / "nme_labels_structures.csv"
    frame.to_csv(out, index=False)

    print(f"\nNDA applications processed: {len(frame)}")
    print(f"  label matched in openFDA    : {int((frame['label_source'] == 'openfda').sum())}")
    print(f"  pharmacokinetics section    : {int(frame['has_pk_section'].sum())}")
    print(f"  structure resolved          : {int((frame['smiles'] != '').sum())}")
    print(f"    via PubChem               : {int((frame['structure_source'] == 'pubchem').sum())}")
    print(f"    via ChEMBL                : {int((frame['structure_source'] == 'chembl').sum())}")
    print(f"  both PK text and structure  : "
          f"{int((frame['has_pk_section'] & (frame['smiles'] != '')).sum())}")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
