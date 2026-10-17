"""Step 5: pull pharmacokinetic parameters out of the label text.

`notes/05` established why this cannot be a plain pattern match. A first-pass
regex for clearance scored 41%, but most hits were *creatinine* clearance from
renal-impairment dosing criteria -- "patients with clearance (CLcr) 15 to 29
mL/min" -- and the rate fell to 16% once that context was excluded. The task is
semantic: separating CL from CL/F from CLcr from renal and dialysis clearance,
and Vss from Vz/F.

So every candidate here is read in the context of its own sentence, each
extraction records the sentence it came from, and anything whose sentence
carries a disqualifying context is rejected with the reason kept rather than
silently dropped. The accompanying validation script measures how often that
judgement is right against a hand-read sample.

Priority follows the feasibility analysis: protein binding and half-life are
the two endpoints with clean units and no apparent-parameter analogue, and are
extracted for use. Clearance and volume are extracted but flagged, because an
apparent value inflated by 1/F cannot be compared with TDC's intravenous
quantities and the prior art went wrong precisely here.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

EXTERNAL_ROOT = Path(__file__).resolve().parents[2] / "data" / "external"

NUMBER = r"\d+(?:[.,]\d+)?"

#: Any of these in a sentence disqualifies a clearance figure from being the
#: systemic clearance of the drug.
CLEARANCE_DISQUALIFIERS = (
    "creatinine", "clcr", "crcl", "cockcroft", "glomerular", "egfr", "gfr",
    "renal impairment", "renal function", "dialysis", "dialyzed", "dialysed",
    "renal clearance", "urinary clearance", "hepatic impairment",
)

#: Marks an apparent (oral, bioavailability-scaled) parameter.
APPARENT = re.compile(r"\b(?:cl|v(?:z|ss|d)?)\s*/\s*f\b|\bapparent\b", re.IGNORECASE)


@dataclass
class Extraction:
    parameter: str
    value: float | None
    unit: str
    is_apparent: bool
    censoring: str  # "", "gt", "lt", "approx"
    value_high: float | None  # populated when the label gives a range
    source_quote: str
    rejected_reason: str = ""


def sentences(text: str) -> list[str]:
    """Split on sentence boundaries, keeping decimals and abbreviations intact."""
    protected = re.sub(r"(\d)\.(\d)", r"\1<DOT>\2", text)
    protected = re.sub(r"\b([A-Z])\.", r"\1<DOT>", protected)
    parts = re.split(r"(?<=[.;])\s+(?=[A-Z(])", protected)
    return [p.replace("<DOT>", ".").strip() for p in parts if p.strip()]


def _to_float(raw: str) -> float | None:
    try:
        return float(raw.replace(",", ""))
    except ValueError:
        return None


def _censoring(sentence: str, span_start: int) -> str:
    """Read the few characters before a number for a censoring marker."""
    prefix = sentence[max(0, span_start - 25) : span_start].lower()
    if re.search(r"(>|greater than|more than|at least)\s*$", prefix):
        return "gt"
    if re.search(r"(<|less than|below|under)\s*$", prefix):
        return "lt"
    if re.search(r"(approximately|about|roughly|~|mean of)\s*$", prefix):
        return "approx"
    return ""


def extract_protein_binding(text: str) -> list[Extraction]:
    """Percent bound to plasma protein, including values reported as fu."""
    found: list[Extraction] = []
    for sentence in sentences(text):
        lowered = sentence.lower()
        if not re.search(r"(protein[- ]?bind|bound to (?:human )?plasma|binding to (?:human )?(?:plasma|serum)|unbound fraction|fraction unbound)", lowered):
            continue
        # Binding to a *target* or a *receptor* is a different quantity.
        if re.search(r"bind(?:s|ing)? to (?:the )?(?:receptor|target|dna|tubulin|enzyme)", lowered):
            continue
        # A sentence whose subject is the metabolite reports the metabolite's
        # binding, not the parent's. Hand-checking caught sofpironium taking
        # "the major sofpironium metabolite (BBI-4010) had plasma protein
        # binding around 2.3-3.7%" as its own value. Sentences that name both,
        # such as "X is 96% and N-desmethyl X is 95% bound", are kept, because
        # the first value quoted is the parent's.
        if re.search(r"\bmetabolite\b", lowered.split("bind")[0] if "bind" in lowered else ""):
            if not re.search(r"\b(?:both|and)\b", lowered.split("metabolite")[0][-40:]):
                continue

        # The sentence filter above has already established that this sentence
        # is about protein binding. Requiring the word "binding" to sit
        # adjacent to the number as well was the first version's mistake: in
        # practice the drug name and a qualifier intervene -- "protein binding
        # of neratinib was greater than 99%" -- and recall fell to a third of
        # what the label text actually contains. So percentages are taken from
        # a binding sentence by default, and only the specific non-binding uses
        # of a percentage are excluded.
        #
        # Label PK sections carry tables flattened into running text, and a
        # table "sentence" can mention binding in one cell and an unrelated
        # percentage in another. Hand-checking 22 extractions found three of
        # five errors arising this way, so a sentence dense in numbers is
        # treated as a table and the percentage must then sit close to the
        # binding word to be believed.
        looks_tabular = len(re.findall(NUMBER, sentence)) > 6 or len(sentence) > 300
        binding_spans = [m.start() for m in re.finditer(r"bind|bound|unbound", lowered)]

        for match in re.finditer(
            rf"({NUMBER})\s*%\s*(?:to|-|–)\s*({NUMBER})\s*%|({NUMBER})\s*(?:to|-|–)\s*({NUMBER})\s*%|({NUMBER})\s*%",
            sentence,
        ):
            groups = match.groups()
            if groups[4] is not None:
                low = high = _to_float(groups[4])
            elif groups[0] is not None:
                low, high = _to_float(groups[0]), _to_float(groups[1])
            else:
                low, high = _to_float(groups[2]), _to_float(groups[3])
            if low is None or high is None or not 0.0 <= low <= 100.0:
                continue

            # A percentage must sit near a binding word, always. Hand-checking
            # a second sample showed the tabular test alone was not enough:
            # long flattened tables (bioequivalence ratios, exposure tables,
            # food-effect statements) still carried a stray percentage far from
            # the binding mention. 120 characters comfortably spans the real
            # phrasings -- "protein binding of neratinib was greater than 99%"
            # puts 37 characters between the two -- while excluding a value
            # from a different cell of the same table.
            reach = 60 if looks_tabular else 120
            if not any(abs(b - match.start()) < reach for b in binding_spans):
                continue

            window = sentence[max(0, match.start() - 55) : match.end() + 45].lower()
            if re.search(
                r"confidence|\bci\b|coefficient of variation|\bcv\b"
                r"|%\s*of the (?:dose|radioactivity|administered)|recovered|excreted"
                r"|increas|decreas|reduc|higher|lower|\bauc\b|c\s*max|bioavailab",
                window,
            ):
                continue

            unbound = bool(
                re.search(r"(?:unbound|free)(?:\s+fraction)?[^.]{0,40}$", sentence[: match.start()].lower())
                or re.search(r"^\s*(?:unbound|free)", sentence[match.end() : match.end() + 20].lower())
            )
            value = (low + high) / 2.0
            found.append(
                Extraction(
                    parameter="ppb",
                    value=(100.0 - value) if unbound else value,
                    unit="% bound",
                    is_apparent=False,
                    censoring=_censoring(sentence, match.start()),
                    value_high=(100.0 - low) if unbound else high,
                    source_quote=sentence[:300],
                )
            )
            break  # one binding value per sentence

        # fu reported as a fraction rather than a percentage.
        for match in re.finditer(rf"(?:fu|fraction unbound|unbound fraction)\s*(?:is|was|of|=|:)?\s*({NUMBER})\b(?!\s*%)", sentence, re.IGNORECASE):
            value = _to_float(match.group(1))
            if value is None or not 0.0 < value < 1.0:
                continue
            found.append(
                Extraction(
                    parameter="ppb",
                    value=100.0 * (1.0 - value),
                    unit="% bound",
                    is_apparent=False,
                    censoring=_censoring(sentence, match.start()),
                    value_high=None,
                    source_quote=sentence[:300],
                )
            )
    return found


def extract_half_life(text: str) -> list[Extraction]:
    """Terminal half-life, normalised to hours."""
    scale = {"hour": 1.0, "hr": 1.0, "h": 1.0, "day": 24.0, "minute": 1 / 60.0, "min": 1 / 60.0}
    found: list[Extraction] = []

    for sentence in sentences(text):
        lowered = sentence.lower()
        if not re.search(r"half[- ]?life|t\s*1\s*/\s*2|t½|terminal elimination", lowered):
            continue
        # Half-life of a metabolite or of an effect is not the drug's t1/2.
        if re.search(r"half[- ]?life of (?:the )?(?:metabolite|m\d|effect)", lowered):
            continue

        pattern = rf"({NUMBER})\s*(?:(?:to|-|–|and)\s*({NUMBER})\s*)?(hours?|hrs?|h|days?|minutes?|min)\b"
        for match in re.finditer(pattern, sentence, re.IGNORECASE):
            low = _to_float(match.group(1))
            high = _to_float(match.group(2)) if match.group(2) else None
            if low is None:
                continue
            unit_word = match.group(3).lower().rstrip("s")
            factor = scale.get(unit_word, scale.get(unit_word[:3]))
            if factor is None:
                continue
            value = ((low + high) / 2.0 if high is not None else low) * factor
            if not 0.01 <= value <= 2000.0:
                continue
            found.append(
                Extraction(
                    parameter="half_life",
                    value=value,
                    unit="h",
                    is_apparent=bool(APPARENT.search(sentence)),
                    censoring=_censoring(sentence, match.start()),
                    value_high=(high * factor) if high is not None else None,
                    source_quote=sentence[:300],
                )
            )
            break  # first numeric half-life in the sentence is the drug's
    return found


def extract_clearance(text: str) -> list[Extraction]:
    """Systemic clearance, with every disqualifying context recorded."""
    found: list[Extraction] = []
    for sentence in sentences(text):
        lowered = sentence.lower()
        if "clearance" not in lowered and not re.search(r"\bcl\s*/?\s*f?\b", lowered):
            continue

        disqualifier = next((d for d in CLEARANCE_DISQUALIFIERS if d in lowered), "")
        pattern = rf"({NUMBER})\s*(l/h(?:r)?(?:/kg)?|ml/min(?:/kg)?|l/hour)\b"
        for match in re.finditer(pattern, sentence, re.IGNORECASE):
            value = _to_float(match.group(1))
            if value is None:
                continue
            found.append(
                Extraction(
                    parameter="clearance",
                    value=value,
                    unit=match.group(2).lower(),
                    is_apparent=bool(APPARENT.search(sentence)),
                    censoring=_censoring(sentence, match.start()),
                    value_high=None,
                    source_quote=sentence[:300],
                    rejected_reason=f"context: {disqualifier}" if disqualifier else "",
                )
            )
    return found


def extract_volume(text: str) -> list[Extraction]:
    """Volume of distribution, flagged when apparent."""
    found: list[Extraction] = []
    for sentence in sentences(text):
        lowered = sentence.lower()
        if not re.search(r"volume of distribution|\bvd\b|\bvss\b|\bvz\b", lowered):
            continue
        for match in re.finditer(rf"({NUMBER})\s*(l/kg|liters?/kg|l\b|liters?)\b", sentence, re.IGNORECASE):
            value = _to_float(match.group(1))
            if value is None:
                continue
            unit = "L/kg" if "kg" in match.group(2).lower() else "L"
            found.append(
                Extraction(
                    parameter="vd",
                    value=value,
                    unit=unit,
                    is_apparent=bool(APPARENT.search(sentence)),
                    censoring=_censoring(sentence, match.start()),
                    value_high=None,
                    source_quote=sentence[:300],
                    rejected_reason="" if unit == "L/kg" else "reported in L, needs a body-weight assumption",
                )
            )
    return found


EXTRACTORS = {
    "ppb": extract_protein_binding,
    "half_life": extract_half_life,
    "clearance": extract_clearance,
    "vd": extract_volume,
}


#: Generic suffixes shared by many INNs, too weak on their own to identify a
#: drug name in a sentence.
_INN_STEMS = re.compile(r"(?:mab|nib|tinib|ciclib|afenib|sartan|prazole|statin|gliflozin|gliptin)\b")


def flag_wrong_drug(extraction: Extraction, ingredient: str) -> Extraction:
    """Mark a value whose sentence names a different drug than this product's.

    Combination products state each component's parameters in turn, and a
    hand-check caught avutometinib taking defactinib's binding value. When the
    sentence names a drug-like word that is not the ingredient, the value is
    kept but flagged rather than discarded, since the sentence may legitimately
    mention a metabolite or a comparator.
    """
    quote = extraction.source_quote.lower()
    stem = ingredient.split()[0].lower() if ingredient else ""
    if not stem or stem in quote:
        return extraction

    others = {w for w in _INN_STEMS.findall(quote)} | {
        w for w in re.findall(r"\b[a-z]{6,}(?:nib|mab|ciclib|afenib)\b", quote)
    }
    if others:
        extraction.rejected_reason = (
            extraction.rejected_reason or "sentence names a different drug"
        )
    return extraction


def extract_all(text: str, ingredient: str = "") -> list[Extraction]:
    if not isinstance(text, str) or not text.strip():
        return []
    results: list[Extraction] = []
    for extractor in EXTRACTORS.values():
        results.extend(extractor(text))
    return [flag_wrong_drug(r, ingredient) for r in results]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    frame = pd.read_csv(EXTERNAL_ROOT / "external_candidates.csv")
    clean = frame[~frame["in_any_tdc_set"] & frame["has_pk_section"]]
    log.info("extracting from %d temporally clean PK sections", len(clean))

    rows = []
    for _, record in clean.iterrows():
        for item in extract_all(record["pk_text"], str(record["primary_ingredient"])):
            rows.append(
                {
                    "application_number": record["application_number"],
                    "ingredient": record["primary_ingredient"],
                    "approval_date": record["approval_date"],
                    "parent_smiles": record["parent_smiles"],
                    "inchikey": record["inchikey"],
                    **asdict(item),
                }
            )

    extractions = pd.DataFrame(rows)
    extractions.to_csv(EXTERNAL_ROOT / "pk_extractions.csv", index=False)

    print(f"\nPK sections processed : {len(clean)}")
    print(f"raw extractions       : {len(extractions)}")
    print("\n=== Compounds with at least one usable value, by parameter ===")
    usable = extractions[extractions["rejected_reason"] == ""]
    for parameter in EXTRACTORS:
        subset = usable[usable["parameter"] == parameter]
        direct = subset[~subset["is_apparent"]]
        print(
            f"  {parameter:10s} any {subset['application_number'].nunique():4d}"
            f"   non-apparent {direct['application_number'].nunique():4d}"
        )

    print("\n=== Clearance rejections, by context ===")
    rejected = extractions[
        (extractions["parameter"] == "clearance") & (extractions["rejected_reason"] != "")
    ]
    print(rejected["rejected_reason"].value_counts().head(10).to_string() or "  none")
    print(f"\nWrote {EXTERNAL_ROOT / 'pk_extractions.csv'}")


if __name__ == "__main__":
    main()
