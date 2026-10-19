"""Replace a section of a markdown draft with the contents of a fragment file.

Written because four attempts to do this through a shell heredoc corrupted the
text: backslashes in LaTeX became control characters, and apostrophes broke the
quoting. Passing prose through a shell is the problem, so the prose lives in
files and only paths cross the boundary.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def splice(target: Path, fragment: Path, start: str, end: str) -> None:
    text = target.read_text(encoding="utf-8")
    if start not in text:
        raise SystemExit(f"start marker not found: {start!r}")
    if end not in text:
        raise SystemExit(f"end marker not found: {end!r}")
    head, _, rest = text.partition(start)
    _, _, tail = rest.partition(end)
    target.write_text(
        head + fragment.read_text(encoding="utf-8").rstrip() + "\n\n" + end + tail,
        encoding="utf-8",
    )
    print(f"spliced {fragment.name} into {target.name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target")
    parser.add_argument("fragment")
    parser.add_argument("--start", required=True, help="first line of the section to replace")
    parser.add_argument("--end", required=True, help="first line of the section that follows")
    args = parser.parse_args()
    splice(Path(args.target), Path(args.fragment), args.start, args.end)


if __name__ == "__main__":
    main()
