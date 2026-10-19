# Thông báo cho TDC trước khi nộp bài A

Ngày soạn: 2026-10-19. **ĐÃ GỬI: issue #405, 2026-10-19**
https://github.com/mims-harvard/TDC/issues/405
Người gửi: MinhNhat-2504. Chưa có phản hồi tính tới lúc ghi.

Lý do phải làm: bài A báo cáo một khiếm khuyết trong công cụ mà cộng đồng đang dùng.
Công bố mà không báo trước cho người bảo trì là sai quy chuẩn, và họ có thể có lý do mình
chưa biết. Ghi lại ở đây để hồ sơ có dấu vết đã thông báo.

---

## Tiêu đề issue

```
Clearance_Hepatocyte_AZ appears to merge the human and rat hepatocyte assays from CHEMBL3301361
```

## Nội dung issue — chép nguyên khối dưới đây

```markdown
Hi, and thanks for maintaining TDC — this report comes from using it seriously
rather than from a drive-by.

While tracing `Clearance_Hepatocyte_AZ` back to its source deposition I found
something I could not reconcile with the dataset description, and a second,
smaller issue alongside it. I may be missing context, so I have written both as
questions with reproducible checks rather than as claims.

## 1. The hepatocyte clearance set appears to contain rat data

The AstraZeneca deposition (ChEMBL `CHEMBL3301361`) holds three relevant
assays, each labelled by species at source:

| Assay | `assay_organism` | n |
|---|---|---|
| `CHEMBL3301370` | *Homo sapiens* — liver microsomes | 1102 |
| `CHEMBL3301372` | *Homo sapiens* — hepatocytes | **408** |
| `CHEMBL3301371` | ***Rattus norvegicus*** — hepatocytes | **837** |

PubChem agrees: AID 1159396 is "Intrinsic clearance measured in **human**
hepatocytes", 408 compounds; AID 1159394 is human liver microsomes, 1102.

`Clearance_Hepatocyte_AZ` has 1213 records over 1020 molecules. Matching each
record against both hepatocyte assays **by value**:

| | records | share |
|---|---|---|
| matches the human hepatocyte value | 370 | 30.5% |
| matches the **rat** hepatocyte value | **796** | **65.6%** |
| both assays report the same value | 31 | 2.6% |
| matches neither | 16 | 1.3% |

`Clearance_Microsome_AZ` (1102) and `PPBR_AZ` (1614) match a single human assay
exactly, so this appears specific to the hepatocyte endpoint.

Reproducible check:

```python
import requests, time
from rdkit import Chem
from tdc.single_pred import ADME

def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else None

def fetch(aid):
    out, off = [], 0
    while True:
        r = requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            params={"assay_chembl_id": aid, "limit": 1000, "offset": off,
                    "only": "canonical_smiles,standard_value"}, timeout=120).json()
        out += r["activities"]
        if len(out) >= r["page_meta"]["total_count"]:
            break
        off += 1000; time.sleep(0.3)
    return {canon(a["canonical_smiles"]): float(a["standard_value"])
            for a in out if a["standard_value"] is not None}

human = fetch("CHEMBL3301372")   # Homo sapiens hepatocytes
rat   = fetch("CHEMBL3301371")   # Rattus norvegicus hepatocytes
tdc   = ADME(name="Clearance_Hepatocyte_AZ").get_data()

h = r = both = neither = 0
for _, row in tdc.iterrows():
    s, v = canon(row["Drug"]), float(row["Y"])
    in_h = s in human and abs(human[s] - v) < 1e-6
    in_r = s in rat   and abs(rat[s]   - v) < 1e-6
    if in_h and in_r: both += 1
    elif in_h: h += 1
    elif in_r: r += 1
    else: neither += 1
print(h, r, both, neither)   # 370 796 31 16
```

The dataset description says the data "contains clearance measures from two
experiments types, hepatocyte and microsomes… we separate them", citing Di et
al. (2012) on human microsomes versus human hepatocytes. That is a rationale
for separating assay *systems*. I could not find anywhere that species is
mentioned.

**Questions:**
- Is the inclusion of `CHEMBL3301371` intentional? If so, could the dataset
  description say so, since models are being trained and ranked on it as a
  human ADME endpoint?
- If not, would you consider either restricting the endpoint to
  `CHEMBL3301372` or adding a species column so users can filter?

I want to be clear that I do not know how this arose and am not assuming it was
careless — pooling by assay type is a reasonable-looking rule if the species
field is not carried through the conversion.

## 2. Measurement qualifiers are dropped in conversion

A smaller and probably related point. The deposition declares its assay ranges
in the assay descriptions themselves:

- `CHEMBL3301372` / `CHEMBL3301370`: "Experimental range **<3 to >150**
  microL/min/1E6 cells" (and per mg for microsomes)
- `CHEMBL3301365`: "Experimental range **10% to 99.95% bound**"

ChEMBL records which values are at a limit in `standard_relation` — 358
microsome and 119 hepatocyte records carry `<` or `>`. In TDC, `Y` is a float
column with no inequality information, and those values appear as exactly 3.0,
150.0 and 99.95.

The result is that 27.4% of `Clearance_Hepatocyte_AZ` and 33.7% of
`Clearance_Microsome_AZ` records sit on a detection limit, indistinguishable
from measurements. A user training a regression model on these has no way to
know.

```python
import numpy as np
from tdc.single_pred import ADME
y = ADME(name="Clearance_Microsome_AZ").get_data()["Y"].to_numpy()
print((y == 3.0).sum(), (y == 150.0).sum(), len(y))   # 287 84 1102
```

**Question:** would you consider carrying `standard_relation` through as an
additional column where the source provides it? It would not change any
existing value or break any existing code, and it would let users choose how to
handle censored data rather than having the choice made silently.

## Context

I found these while preparing a manuscript on the pharmacokinetic endpoints,
and I am raising them here before submitting rather than after. If either turns
out to be intentional or already documented somewhere I have missed, I would
rather correct the manuscript than publish something misleading — so please do
say.

Happy to open a PR for a species column or a qualifier column if that would be
useful, and happy to share the full analysis and scripts.
```

---

## Hướng dẫn gửi

### Bước 1 — Mở issue mới

Vào https://github.com/mims-harvard/TDC/issues và bấm **New issue**.
Nếu nó hiện mẫu chọn loại issue, chọn **Blank issue** (hoặc "Open a blank issue"
ở dòng chữ nhỏ bên dưới) — đây không phải báo lỗi phần mềm nên các mẫu bug report
không hợp.

### Bước 2 — Trước khi dán, tìm xem đã có ai nêu chưa

Tìm trong issue của repo: `hepatocyte`, `rat`, `species`, `CHEMBL3301371`, `clearance`.
Nếu đã có issue tương tự thì **bình luận vào đó** thay vì mở cái mới — vừa lịch sự vừa
dễ được trả lời hơn.

### Bước 3 — Dán

Tiêu đề và nội dung lấy từ hai khối trên. **Khối nội dung đã là markdown**, dán thẳng
vào ô body, GitHub sẽ render bảng và code block đúng.

Lưu ý khi dán: khối code Python nằm **lồng trong** khối markdown ngoài. Khi chép, chép
phần **bên trong** cặp ```` ```markdown ```` và ```` ``` ```` ngoài cùng — đừng chép cả
hai dấu ngoặc đó.

### Bước 4 — Sau khi gửi

- Dán link issue vào đây để mình ghi vào bài (Mục 5.5 hiện ghi "đã thông báo", cần số
  issue để trích dẫn được).
- **Chờ một khoảng hợp lý trước khi nộp bài.** Thông lệ là 2–4 tuần. Nếu họ trả lời
  nhanh thì tốt; nếu im lặng thì vẫn nộp được, nhưng trong bài phải ghi rõ: *"reported to
  the maintainers as issue #N on <ngày>; no response at time of submission."*
- Nếu họ trả lời rằng việc gộp là **có chủ đích**, bài phải sửa: luận điểm chuyển từ
  "không khai báo" sang "khai báo chưa đủ rõ", vẫn đứng được nhưng nhẹ hơn.
- Nếu họ **sửa dataset**, đó là kết quả tốt nhất — và bài vẫn có giá trị vì nó ghi lại
  trạng thái tại thời điểm công bố, kèm số phiên bản.

### Bước 5 — Ghi vào bài

Sau khi có số issue, Mục 5.5 sẽ ghi:

> We reported both observations to the TDC maintainers (issue #N, YYYY-MM-DD)
> before submission. `[trạng thái phản hồi]`

---

## Hai điều cân nhắc trước khi bấm gửi

1. **Issue công khai.** Ai cũng đọc được, kể cả trước khi bài đăng. Nếu bạn lo bị lấy ý
   tưởng thì lưu ý nội dung issue chỉ nêu **hai quan sát về dữ liệu**, không nêu phần
   phân tích hệ quả (trần nhiễu bị hỏng, fu sai 8,3 lần) — đó mới là phần khó và là đóng
   góp của bài.
2. **Giọng điệu.** Bản nháp cố ý nhẹ nhàng và đặt câu hỏi. Đừng sửa thành khẳng định
   mạnh hơn: nếu hoá ra họ có lý do chính đáng thì giọng khẳng định sẽ phản tác dụng,
   còn nếu không có thì giọng nhẹ vẫn đủ sức nặng vì số liệu tự nói.
