"""Renumber references by order of first citation, computing the order.

ACS numbers references in order of first appearance. Keeping that correct by
hand failed twice while this paper was being assembled -- once when a removed
entry shifted everything below it, and once when two sequential replacements
collided and sent two citations to the same marker.

So the numbering is derived rather than maintained. Each reference has a stable
key; the text may carry either `[[key]]` or a rendered superscript, and this
reads both, assigns numbers in order of first appearance, and rewrites the
markers and the list together.

It is idempotent: the keyed block is written in the same order as the rendered
list, so a rendered superscript can always be mapped back to its key on the
next run. The first version was not, which is how it broke the first time a
reference was added after rendering.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper" / "paper_A_benchmark_defects.md"

DIGITS = "⁰¹²³⁴⁵⁶⁷⁸⁹"
SUPERSCRIPT = {str(d): c for d, c in enumerate(DIGITS)}
PLAIN = {c: str(d) for d, c in enumerate(DIGITS)}
TOKEN = re.compile(r"\[\[([a-z0-9_]+)\]\]")
#: A run of adjacent tokens is one citation group and renders separated by
#: an apostrophe; joining the digits directly would make 4 and 5 read as 45.
TOKEN_RUN = re.compile(r"(?:\[\[[a-z0-9_]+\]\])+")
RUN = re.compile(rf"[{DIGITS}]+(?:'[{DIGITS}]+)*")
MARKER = "<!-- REFERENCES -->"


def to_superscript(number: int) -> str:
    return "".join(SUPERSCRIPT[d] for d in str(number))


def load_entries(block: str) -> dict[str, str]:
    """Keyed reference bodies, in file order."""
    entries: dict[str, str] = {}
    key, buffer = None, []
    for line in block.split("\n"):
        match = re.match(r"^@([a-z0-9_]+):\s*(.*)$", line)
        if match:
            if key:
                entries[key] = " ".join(buffer)
            key, buffer = match.group(1), [match.group(2)]
        elif key is not None:
            buffer.append(line.strip())
    if key:
        entries[key] = " ".join(buffer)
    return {k: re.sub(r"\s+", " ", v).strip() for k, v in entries.items()}


def main() -> None:
    text = PAPER.read_text(encoding="utf-8")
    head, _, block = text.partition(MARKER)
    entries = load_entries(block)
    keys = list(entries)

    # Drop any previously rendered list; it is regenerated below.
    body = head.split("## References")[0]

    # Rendered superscripts map back to keys by position, which holds because
    # this script always writes the keyed block in rendered order.
    def detokenise(match: re.Match) -> str:
        out = []
        for part in match.group(0).split("'"):
            index = int("".join(PLAIN[c] for c in part)) - 1
            if not 0 <= index < len(keys):
                raise SystemExit(f"superscript {index + 1} has no reference")
            out.append(f"[[{keys[index]}]]")
        return "".join(out)

    body = RUN.sub(detokenise, body)

    order: list[str] = []
    for key in TOKEN.findall(body):
        if key not in entries:
            raise SystemExit(f"cited but not defined: [[{key}]]")
        if key not in order:
            order.append(key)

    unused = [k for k in entries if k not in order]
    if unused:
        raise SystemExit(f"defined but never cited: {', '.join(unused)}")

    numbers = {key: i + 1 for i, key in enumerate(order)}
    def render(match: re.Match) -> str:
        found = TOKEN.findall(match.group(0))
        return "'".join(to_superscript(numbers[k]) for k in found)

    body = TOKEN_RUN.sub(render, body)

    rendered = "## References\n\n" + "\n\n".join(
        f"({numbers[key]}) {entries[key]}" for key in order
    )
    keyed = MARKER + "\n\n" + "\n\n".join(f"@{key}: {entries[key]}" for key in order)

    PAPER.write_text(
        body.rstrip() + "\n\n" + rendered + "\n\n" + keyed + "\n", encoding="utf-8"
    )

    print(f"{len(order)} references, numbered by first appearance:")
    for key in order:
        print(f"  {numbers[key]:2d}  {key}")


if __name__ == "__main__":
    main()
