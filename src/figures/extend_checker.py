"""Let `check_paper.py` check either manuscript.

The checker was written for paper A and names it as a constant. Paper C has the
same failure modes -- a stale cross-reference, a figure referenced but not
rendered, a citation defined and never used -- and nothing was checking them.

The number check stays specific to paper A, since the figures it verifies are
paper A's and the supporting information belongs to it. The three structural
checks apply to whichever file is named on the command line.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "src" / "figures" / "check_paper.py"

IMPORTS = ("""from __future__ import annotations

import re
from pathlib import Path""",
"""from __future__ import annotations

import re
import sys
from pathlib import Path""")

CONSTANT = ("""PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md\"""",
            """DEFAULT_PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md\"""")

MAIN = ('''def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    si = SI.read_text(encoding="utf-8")

    sections = {
        "numbers against source files": check_numbers(text + si),
        "cross-references": check_cross_references(text),
        "figures": check_figures(text),
        "citations": check_citations(text),
        "retracted figures": check_retractions(text),
    }
''',
'''def main() -> None:
    paper = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PAPER
    if not paper.is_absolute():
        paper = ROOT / paper
    text = paper.read_text(encoding="utf-8")
    is_paper_a = paper.name == DEFAULT_PAPER.name

    sections = {
        "cross-references": check_cross_references(text),
        "figures": check_figures(text),
        "citations": check_citations(text),
    }
    if is_paper_a:
        # The quoted figures and the retraction passages are paper A's, and the
        # supporting information belongs to it.
        si = SI.read_text(encoding="utf-8")
        sections = {
            "numbers against source files": check_numbers(text + si),
            **sections,
            "retracted figures": check_retractions(text),
        }
    print(f"checking {paper.name}")
''')


def main() -> None:
    text = CHECKER.read_text(encoding="utf-8")
    if "DEFAULT_PAPER" in text:
        print("already parameterised")
        return
    for old, new in (IMPORTS, CONSTANT, MAIN):
        if text.count(old) != 1:
            raise SystemExit(f"anchor not unique or missing:\n  {old[:70]}")
        text = text.replace(old, new, 1)
    CHECKER.write_text(text, encoding="utf-8")
    print("check_paper.py now takes a paper path")


if __name__ == "__main__":
    main()
