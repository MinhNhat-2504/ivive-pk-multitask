# Finding 15 — `Clearance_Hepatocyte_AZ` của TDC gộp dữ liệu **người và chuột** vào một endpoint

Ngày: 2026-10-17. Truy vấn trực tiếp ChEMBL API, đối chiếu từng giá trị với TDC.

> **Phát hiện này đồng thời (a) là kết quả mạnh nhất của dự án và (b) phá huỷ phân tích
> trần nhiễu ở `notes/03` và `notes/04`.** Cả hai đều phải viết lại.

## 1. Bằng chứng

Bản nộp AstraZeneca `CHEMBL3301361` có 15 assay. Ba cái liên quan:

| Assay | Loài | n |
|---|---|---|
| `CHEMBL3301370` | **Homo sapiens** — microsome gan | 1 102 |
| `CHEMBL3301372` | **Homo sapiens** — tế bào gan | **408** |
| `CHEMBL3301371` | **Rattus norvegicus** — tế bào gan | **837** |

TDC `Clearance_Hepatocyte_AZ` có **1 213 bản ghi / 1 020 phân tử**. Còn 408 + 837 = **1 245**.

### Khớp theo giá trị (không phải chỉ theo phân tử), trên 1 213 bản ghi TDC

| | n | % |
|---|---|---|
| Khớp giá trị assay **NGƯỜI** | 370 | 30,5% |
| **Khớp giá trị assay CHUỘT** | **796** | **65,6%** |
| Hai assay trùng giá trị | 31 | 2,6% |
| Không khớp | 16 | 1,3% |

> **Khoảng hai phần ba `Clearance_Hepatocyte_AZ` là dữ liệu tế bào gan CHUỘT**, nằm trong
> một benchmark ADMET được trình bày như endpoint trên người.

Hai endpoint kia **sạch** — khớp chính xác một assay người duy nhất:
`cl_mic` = 1 102 = `CHEMBL3301370`; `ppbr` = 1 614 = `CHEMBL3301365`.

## 1b. Việc gộp xảy ra ở TDC, không phải thượng nguồn — đã truy đủ

**Nguồn thượng nguồn sạch và ghi đúng loài:**

| Nguồn | Mô tả | n |
|---|---|---|
| PubChem **AID 1159396** | *"Intrinsic clearance measured in **human** hepatocytes"* | **408** |
| PubChem AID 1159394 | *"...**human** liver microsomes"* | 1 102 |
| ChEMBL `CHEMBL3301372` | `assay_organism: Homo sapiens` | 408 |
| ChEMBL `CHEMBL3301371` | `assay_organism: **Rattus norvegicus**` | 837 |

**Tài liệu của TDC không khai báo việc gộp loài.** Nguyên văn từ tdcommons.ai:

> *"It contains clearance measures from **two experiments types, hepatocyte and
> microsomes**. As many studies have shown various clearance outcomes given these two
> different types, **we separate them**."*
> *"Dataset Statistics: 1,102/1,020 drugs for microsome/hepatocyte clearance."*

Tách theo **loại thí nghiệm**, **không nhắc tới loài**. Và trích dẫn biện minh là Di et
al. 2012, *"comparing intrinsic clearance values between **human** liver microsomes and
hepatocytes"* — một bài so người với người, dùng để biện minh cho một bộ dữ liệu mà hai
phần ba là chuột.

→ Quy trách nhiệm được: **ChEMBL và PubChem ghi đúng; TDC trộn hai assay khi curate và
mô tả phép tách chỉ theo loại thí nghiệm.** Mọi mô hình huấn luyện trên TDC kế thừa điều
này mà không có cách nào biết.

## 2. ⚠ Hệ quả: trần nhiễu ở `notes/03` và `notes/04` đo nhầm đại lượng

Ở `notes/03` mình tìm thấy **193 "nhóm đo lặp"** trong `cl_hep` và dùng độ tản trong nhóm
để ước lượng độ lặp lại của phép đo.

**187 trong 193 nhóm đó (97%) đúng bằng một cặp giá trị người + chuột của cùng một phân tử.**

Nghĩa là `sd_within` = 1,024 log (hệ số 2,8 lần) **không phải nhiễu phép đo**. Nó là
**khác biệt loài giữa tế bào gan người và chuột** — một đại lượng sinh học có thật, đã
biết từ lâu, và hoàn toàn không phải thứ mình tưởng đang đo.

Vậy:

| Con số | Trạng thái |
|---|---|
| ICC 0,294, trần 0,541 (`notes/03`) | ❌ **rút lại** — đo tương hợp giữa loài |
| ICC 0,571, trần ≥0,764 (`notes/04`) | ❌ **rút lại** — cùng lỗi, chỉ bỏ bớt kiểm duyệt |
| Kết luận "`cl_hep` chưa bão hoà, còn dư địa ≥0,22" | ❌ **không có cơ sở** |

**Không còn ước lượng trần nhiễu nào cho `cl_hep`.** Và không dựng lại được từ dữ liệu
công khai: assay người có 408 bản ghi trên 408 phân tử — **không đo lặp nào**.

Sai lầm cũng giải thích luôn một chi tiết lúc đó mình thấy lạ mà không truy đến cùng:
**0 nhóm nào trùng khít giá trị**. Người và chuột hiếm khi cho đúng cùng một số. Mình đã
diễn giải nó thành "TDC gộp bỏ các cặp trùng khớp" — sai.

## 3. Phần nào ở `notes/04` vẫn đứng

**Kiểm duyệt vẫn thật và vẫn đo đúng.** Nó là tính chất của phân bố giá trị, không phụ
thuộc vào việc bản ghi đến từ loài nào:

- `cl_hep`: 195 bản ghi tại đúng 3,0 và 137 tại đúng 150,0 → **27,4%**
- `cl_mic`: 287 và 84 → **33,7%**

Và giờ có nguồn trích dẫn trực tiếp, không cần suy từ histogram nữa. Mô tả assay của chính
bản nộp ghi nguyên văn:

> `CHEMBL3301372`: *"Intrinsic clearance measured in human hepatocytes following incubation
> at 37C. **Experimental range <3 to >150 microL/min/1E6 cells.**"*
> `CHEMBL3301370`: *"...**Experimental range <3 to >150 microL/min/mg.**"*
> `CHEMBL3301365`: *"% bound to plasma by equilibrium dialysis... **Experimental range 10%
> to 99.95% bound.**"*

**ChEMBL giữ cờ kiểm duyệt trong trường `standard_relation`; TDC làm mất.** Microsome:
274 `<` + 84 `>` = 32,5%. Thông tin mất **trong khâu chuyển đổi của TDC**, không phải ở
bản nộp gốc.

Thêm: 13 bản ghi ở đúng 3,0 trong `CHEMBL3301370` **không** mang cờ `<`. Nên lọc theo cờ
bỏ sót ~5%; **phát hiện theo giá trị — cách mình làm — bắt được hết.**

## 4. PPBR cũng bị kiểm duyệt, và chưa ai nói — mới

Mình chưa từng xét kiểm duyệt ở `ppbr`. Đo lại:

- max = **đúng 99,95%** = trần assay khai báo. **4 bản ghi ghim ở đó.**
- **299 bản ghi (18,5%) có PPB ≥ 99%** → fu ≤ 0,01
- 187 (11,6%) ≥ 99,5% → fu ≤ 0,005
- **0 bản ghi nào mang cờ `standard_relation`** — kiểm duyệt *về bản chất* nhưng không có
  trong metadata, nên lọc theo cờ sẽ bỏ sót hoàn toàn.

**Điều này quan trọng với mô hình của mình** vì fu = 1 − PPB/100 là biến trung gian của
ràng buộc IVIVE. Sai số trong fu bị khuếch đại theo 1/(1 − PPB/100): ở 99,95% hệ số
khuếch đại là **~2 000 lần**. PPB 99,95% so với 99,90% chỉ lệch 0,05 điểm phần trăm nhưng
là **chênh 2 lần trong fu**.

→ Phải kiểm duyệt hoặc loại PPB ≥ 99,9% thay vì truyền fu qua đó.

## 5. Những việc phải làm lại

- [ ] **Rút lại mọi con số trần nhiễu** khỏi bài (Mục 5.2, 5.4). Thay bằng: không ước
      lượng được từ dữ liệu công khai, và nêu lý do.
- [ ] Viết lại Mục 5 quanh **gộp loài** + **kiểm duyệt**, hai thứ đo được và có nguồn.
- [ ] Kiểm xem việc gộp loài có ảnh hưởng mọi kết quả mô hình không — baseline, Tobit,
      giảm nhãn, OOD đều huấn luyện trên tập lẫn hai loài này.
- [ ] Xử lý kiểm duyệt trần của `ppbr` trong ràng buộc IVIVE.
- [x] ~~Kiểm tra TDC có ghi chú việc gộp này không~~ → **không có**, xem §1b.

## 6. Hệ quả cho toàn bộ dự án — nói thẳng

Mô hình ràng buộc của mình ghép `cl_hep` với `cl_mic` qua một phép hiệu chỉnh affine học
được, và hoá ra `cl_hep` là **hỗn hợp hai loài**. Hệ số hiệu chỉnh Deming học được
(dốc 0,91–1,11) thực chất đang bắc cầu giữa *microsome người* và *hỗn hợp tế bào gan
người-chuột*.

Điều đó **không làm sai** các kết quả so sánh ghép cặp — mọi arm đều thấy cùng dữ liệu,
nên kết luận "học phép hiệu chỉnh thắng cố định nó" vẫn đứng. Nhưng **diễn giải dược học
thay đổi**: phép hiệu chỉnh đang hấp thụ một phần khác biệt loài, không chỉ khác biệt hệ
thống in vitro.

Và nó mở ra một thí nghiệm sắc hơn: **tách 370 bản ghi người ra, huấn luyện riêng, so
sánh.** Nếu ràng buộc vẫn giúp trên tập người thuần thì luận điểm mạnh hơn hẳn.
