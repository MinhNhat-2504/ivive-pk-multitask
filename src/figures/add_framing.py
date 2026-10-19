"""Add two references and the framing one of them supplies.

Nael et al. published a Viewpoint in JCIM in February 2026 arguing that pooling
non-equivalent assays under one label is itself a defect, and calling for
benchmarks to declare their aggregation logic. That is the general claim this
paper gives a concrete case of, on the most-used ADMET benchmark, so it belongs
in the introduction rather than buried in related work.

Williams and Richard make the database-integrity version of the same point and
belong in the related-work paragraph on source consistency.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"

NEW_REFERENCES = {
    "nael": (
        "Nael, M. A.; Alakonda, L. M.; Elokely, K. M. Defining the Data Set Defines "
        "the QSAR Claim. *J. Chem. Inf. Model.* **2026**, *66*, 2951–2954. "
        "DOI: 10.1021/acs.jcim.6c00514."
    ),
    "williams": (
        "Williams, A. J.; Richard, A. M. Three Pillars for Ensuring Public Access and "
        "Integrity of Chemical Databases Powering Cheminformatics. *J. Cheminform.* "
        "**2025**, *17*, 40. DOI: 10.1186/s13321-025-00983-9."
    ),
}

EDITS = [
    # Introduction: the general claim this paper is a case of.
    (
        "duplication. Neither examined what the labels are.",
        "duplication. Neither examined what the labels are.\n\n"
        "A Viewpoint in this journal has since argued that the problem begins "
        "earlier than either: that public resources combine assays of different "
        "design under one label, that doing so \"can merge phenomena that are not "
        "equivalent\", and that a benchmark should therefore declare its endpoint "
        "definitions and aggregation logic explicitly.[[nael]] That is a general "
        "claim about how benchmarks are built. This paper is a concrete case of it, "
        "on the most used ADMET collection, found by going back to the source.",
    ),
    # Related work: the database-integrity statement of the same concern.
    (
        "**Source consistency.** Parrondo-Pizarro et al.",
        "**Database integrity and source consistency.** Williams and Richard argue "
        "that chemical databases are widely trusted without independent validation, "
        "and that aggregating content across sources \"makes it difficult to trace "
        "the origin of data\".[[williams]] The two defects reported here are "
        "instances of exactly that, in a resource built on top of such a "
        "database.\n\nParrondo-Pizarro et al.",
    ),
]


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    for old, new in EDITS:
        if old not in text:
            raise SystemExit(f"anchor not found:\n  {old[:80]}")
        text = text.replace(old, new, 1)

    marker = "<!-- REFERENCES -->\n\n"
    head, _, block = text.partition(marker)
    additions = "\n\n".join(f"@{key}: {body}" for key, body in NEW_REFERENCES.items())
    text = head + marker + block.rstrip() + "\n\n" + additions + "\n"

    PAPER.write_text(text, encoding="utf-8")
    print(f"added {len(NEW_REFERENCES)} references and their framing")


if __name__ == "__main__":
    main()
