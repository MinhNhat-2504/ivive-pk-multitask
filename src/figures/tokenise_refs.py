"""One-off: convert the numbered reference list to the keyed form.

Turns each `(n) Author...` entry into `@key: Author...` in a block at the end
of the file, and each superscript marker into `[[key]]`, so that
`renumber_refs.py` can take over and the numbering stops being something a
human has to track.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"

#: current number -> stable key
KEYS = {
    1: "tdc_neurips", 2: "tdc_natchembio", 3: "koleiev", 4: "schuh",
    5: "parrondo", 6: "chembl_az", 7: "temesi", 8: "di", 9: "obach",
    10: "lombardo", 11: "pksmart", 12: "bowman", 13: "wilson", 14: "gidal",
    15: "tobin", 16: "lightgbm", 17: "rdkit",
}

DIGIT = {c: str(d) for d, c in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹")}
RUN = re.compile(r"[⁰¹²³⁴-⁹]+(?:'[⁰¹²³⁴-⁹]+)*")


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    body, _, reflist = text.partition("## References")

    entries: dict[int, str] = {}
    number = None
    buffer: list[str] = []
    for line in reflist.split("\n"):
        match = re.match(r"^\((\d+)\) (.*)$", line)
        if match:
            if number:
                entries[number] = " ".join(buffer)
            number, buffer = int(match.group(1)), [match.group(2)]
        elif number and line.strip() and not line.startswith("#"):
            buffer.append(line.strip())
        elif number and not line.strip() and buffer:
            entries[number] = " ".join(buffer)
            number, buffer = None, []
    if number and buffer:
        entries[number] = " ".join(buffer)

    missing = set(KEYS) - set(entries)
    if missing:
        raise SystemExit(f"no entry found for reference(s) {sorted(missing)}")

    def tokenise(match: re.Match) -> str:
        out = []
        for part in match.group(0).split("'"):
            value = int("".join(DIGIT[c] for c in part))
            out.append(f"[[{KEYS[value]}]]")
        return "".join(out)

    body = RUN.sub(tokenise, body)

    def normalise(value: str) -> str:
        return re.sub(r"\s+", " ", value).strip()

    block = "<!-- REFERENCES -->\n\n" + "\n\n".join(
        f"@{KEYS[n]}: {normalise(entries[n])}" for n in sorted(KEYS)
    ) + "\n"
    PAPER.write_text(body.rstrip() + "\n\n## References\n\n" + block, encoding="utf-8")
    print(f"tokenised {len(KEYS)} references")


if __name__ == "__main__":
    main()
