"""Fix what a public repository cited by two papers should not have wrong.

The repository is private. It has to be public, because the issue comment and
both manuscripts' code-availability statements point at it, and the Zenodo
record declares it as a supplement. Before flipping it, three things:

  - a link to `paper/draft.md`, which is gitignored and therefore does not
    exist on GitHub at all;
  - a summary line still asserting a corrected noise ceiling of 0.764, which
    paper A retracts outright: the apparent replicate groups turned out to be
    human/rat pairs, so both the original estimate and the correction were
    withdrawn. A repository cited by that paper should not contradict it;
  - no LICENSE, which means default copyright applies and nobody may legally
    reuse code the papers describe as available.

MIT for the code, because that is the convention for research code accompanying
a paper and imposes nothing on anyone. The data files keep their own terms and
the README now says so: the FDA label text is a US government work and the
ChEMBL deposition is CC BY-SA 3.0.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
LICENSE = ROOT / "LICENSE"

MIT = """MIT License

Copyright (c) 2026 Nhat Trinh Ngoc Minh and Linh Nguyen Thi My

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

LICENCE_SECTION = """
## Giấy phép

Mã nguồn: MIT, xem [LICENSE](LICENSE).

Dữ liệu trong `data/` giữ giấy phép riêng của nguồn. Văn bản nhãn thuốc FDA là
tác phẩm của chính phủ Hoa Kỳ; bản gửi AstraZeneca trên ChEMBL (document
`CHEMBL3301361`) là CC BY-SA 3.0. Bộ dữ liệu đã kiểm chứng và phần chú giải
xuất xứ được lưu trữ riêng tại Zenodo dưới CC BY 4.0, concept DOI
[10.5281/zenodo.23462435](https://doi.org/10.5281/zenodo.23462435).

## Bản thảo

Các bản thảo không nằm trong repo này (thư mục `paper/` bị `.gitignore` loại
trừ). Dữ liệu và mã để tái lập mọi con số trong đó thì có.
"""

EDITS = [
    (
        "Xem `notes/` theo thứ tự, và bản thảo paper ở [paper/draft.md](paper/draft.md).",
        "Xem `notes/` theo thứ tự. Bản thảo không nằm trong repo; xem mục Giấy phép\nở cuối.",
    ),
    (
        """- [03 — trần nhiễu](notes/03-noise-ceiling.md): phương pháp ước lượng trần nhiễu.
  **⚠ Con số trong trang này đã bị `notes/04` bác bỏ một phần.**
- [04 — kiểm duyệt assay](notes/04-censoring.md): 27–34% nhãn clearance là cờ ngưỡng phát
  hiện, không phải phép đo. Sửa lại trần nhiễu từ 0,541 lên **≥ 0,764** — tức endpoint
  **chưa bão hoà**, vẫn còn dư địa thật.""",
        """- [03 — trần nhiễu](notes/03-noise-ceiling.md): phương pháp ước lượng trần nhiễu.
  **⚠ Toàn bộ hướng này đã bị thu hồi.** Xem `notes/12`: các "nhóm lặp" dùng để
  ước lượng hoá ra là cặp người/chuột, nên cả con số gốc 0,541 lẫn con số "đã
  sửa" 0,764 đều không có cơ sở. Giữ lại hai trang này vì bản thảo mô tả chính
  lần sai đó.
- [04 — kiểm duyệt assay](notes/04-censoring.md): 27–34% nhãn clearance là cờ ngưỡng phát
  hiện, không phải phép đo. Phần kiểm duyệt đứng vững; phần hiệu chỉnh trần nhiễu
  trong trang này thì không, xem trên.""",
    ),
]


def main() -> None:
    if not LICENSE.exists():
        LICENSE.write_text(MIT, encoding="utf-8")
        print("LICENSE written (MIT)")
    else:
        print("LICENSE already present, left alone")

    text = README.read_text(encoding="utf-8")
    applied = 0
    for old, new in EDITS:
        if new in text and old not in text:
            continue
        if text.count(old) != 1:
            raise SystemExit(f"README anchor not unique or missing:\n  {old[:80]}")
        text = text.replace(old, new, 1)
        applied += 1

    if "## Giấy phép" not in text:
        text = text.rstrip() + "\n" + LICENCE_SECTION
        applied += 1

    README.write_text(text, encoding="utf-8")
    print(f"README: {applied} change(s)")


if __name__ == "__main__":
    main()
