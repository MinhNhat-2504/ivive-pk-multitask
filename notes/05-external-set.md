# Finding 05 — External set theo thời gian từ FDA: làm được, nhưng không cho nhánh in vitro

Ngày tra: 2026-10-17. Tra cứu nguồn dữ liệu bằng truy vấn trực tiếp. **Chưa dựng dataset**
— trang này là kết quả khảo sát khả thi + kế hoạch.

> Các số dưới đây đến từ truy vấn thật vào openFDA / DailyMed / PubChem / ChEMBL và từ
> việc tải về hai hồ sơ thẩm định đầy đủ. Những chỗ chỉ là suy đoán đều được đánh dấu.

## 1. Kết luận ngắn

| Endpoint TDC | n external ước tính | Kết luận |
|---|---|---|
| `ppbr` (% gắn) | **150–200** | ✅ Validation chính |
| `half_life` (giờ) | **140–165** | ✅ Validation chính |
| `vdss` (L/kg) | 35–55 | ⚠ Phụ, chỉ phần IV thật |
| CL toàn thân in vivo (đầu ra IVIVE) | 30–40 | ⚠ Kiểm tra chuỗi đầu-cuối |
| `cl_hep`, `cl_mic` | **0** | ❌ **Không lấy được** |

## 2. Nhánh in vitro không validate được — đã kiểm chứng, không phải suy đoán

Tải nguyên hai hồ sơ Clinical Pharmacology và chạy trích xuất văn bản:

- **Istradefylline ClinPharmR (2019)**, 362 KB text: `"intrinsic clearance"` → **0**,
  `CLint` → **0**, `µL/min` → **0**. Các lần nhắc tế bào gan/microsome đều về định danh
  CYP, hằng số ức chế tương tác thuốc, IC50 — không phải độ ổn định chuyển hoá.
- **Zuranolone IntegratedR (2023)**, 730 KB text: `"intrinsic clearance"` → **0**,
  `CLint` → **0**.

Zuranolone *có* mục "Stability in Hepatic Microsomes and Hepatocytes", nhưng nó báo
**"predicted hepatocyte clearance ~0,88 L/h/kg"** — tức giá trị **đã quy đổi sang toàn
thân**, không phải số đọc của assay. Muốn nghịch đảo về µL/min/10⁶ tế bào thì cần mật độ
tế bào gan, khối lượng gan, mô hình quy đổi (well-stirred hay parallel-tube), có áp hiệu
chỉnh fu,inc và fu,p hay không, và cơ sở máu hay huyết tương — **hồ sơ không nêu cái nào**.
Phép nghịch đảo đó không xác định được, và sai số của giả định sẽ lớn hơn tín hiệu cần đo.

Cộng với kiểm duyệt ở `notes/04`, nhánh in vitro kháng lại kiểm định ngoài vì **hai lý do
độc lập**. Viết thẳng điều này vào bài — nó làm bài mạnh lên, không yếu đi.

## 3. Vấn đề `/F` — và một lỗi trong công trình gần nhất

Nhãn thuốc báo tham số **biểu kiến** nhiều hơn tham số thật. Trên 264 nhãn:

| | số nhãn |
|---|---|
| Vd biểu kiến (`Vz/F`, `Vss/F`) | 101 |
| Vd không biểu kiến | 59 |
| Vd báo bằng **L/kg** | chỉ 10 (115 báo bằng **L**) |
| `CL/F` | 84 |
| CL toàn thân thật | ~30–35 |

Giá trị biểu kiến bị thổi lên theo 1/F. Với thuốc phân tử nhỏ đường uống hiện đại, F
thường 0,1–0,5 → **thổi 2–10 lần**, theo từng chất, **không hiệu chỉnh được**. Trong khi
`VDss_Lombardo` là VDss tĩnh mạch thật.

**PKSmart (Seal et al., J Cheminform 2025) đã mắc đúng lỗi này.** Bộ external FDA của họ
ghi mobocertinib VDss = **58,48 L/kg** (≈4 100 L) và infigratinib 26,67 L/kg — đó là
`Vz/F` chia 70 kg rồi dán nhãn VDss. VDss 58 L/kg không hợp lý về sinh lý.

→ **Một bộ external theo thời gian được sàng `/F` đúng cách hiện chưa tồn tại.** Dựng
đúng là một đóng góp phương pháp, không phải việc lắp ráp.

Một cảnh báo kỹ thuật: regex cho "clearance" cho tỷ lệ 41%, nhưng phần lớn là **độ thanh
thải creatinin** trong tiêu chí chỉnh liều suy thận. Loại ngữ cảnh đó ra thì còn 16%.
Phải trích xuất **theo ngữ nghĩa**, bắt mô hình xuất kèm: tham số, đơn vị, cờ biểu kiến,
đường xác định (IV/uống/popPK), tình trạng kiểm duyệt, và câu trích nguồn.

## 4. Hai endpoint sạch

- **`ppbr`** — tốt nhất: không thứ nguyên, không cần giả định cân nặng, không có khái
  niệm `/F`, phủ 78% nhãn. Lưu ý: nhiều giá trị **bị kiểm duyệt** (">99,5%", "<10%") nên
  phải xử lý như dữ liệu kiểm duyệt; một số nhãn báo **fu** nên phải đổi %gắn = 100×(1−fu);
  PPBR_AZ là assay nội bộ AstraZeneca còn nhãn lấy từ nhà tài trợ → **chờ có lệch hệ
  thống giữa hai assay**, nên báo cáo tương quan hạng song song với MAE.
- **`half_life`** — đơn vị trực tiếp, phủ 61%. Lưu ý ngày (11 ca) và phút (4 ca), và
  khoảng ("19,7 đến 24,6 giờ" → lấy trung điểm). Obach là t½ tĩnh mạch còn nhãn là đường
  uống; dưới động học flip-flop thì t½ uống phản ánh hấp thu — vấn đề thật nhưng nhỏ hơn
  nhiều so với `/F`.

## 5. Mốc thời gian an toàn

Dữ liệu AZ: ChEMBL **CHEMBL3301361**, nộp ~2015–2016. Lombardo 2016. Obach 2008.
→ **Mốc 2017-01-01 sạch cho cả năm endpoint.** Không cần đẩy mốc muộn hơn.

## 6. Đường ống dựng dữ liệu

1. Tải `drug-drugsfda-0001-of-0001.json.zip` (8,8 MB). Lọc `ORIG` + `AP` +
   `Type 1 - New Molecular Entity`, ngày ≥ 2017-01-01 → **451**. Giữ `NDA*` → **324**.
2. Lấy nhãn: openFDA `/drug/label.json`, dùng trường **`pharmacokinetics`** (SPL 12.3),
   *không* phải `clinical_pharmacology`. ~19% trượt → dự phòng DailyMed
   `spls.json?drug_name=` rồi SPL XML, nhắm **LOINC 43682-4** (không phải 34090-1 — đó
   là mục cha).
3. SMILES: PubChem PUG-REST theo tên, dùng **`IsomericSMILES`** (`CanonicalSMILES` của
   PubChem **bỏ mất lập thể**); dự phòng ChEMBL. Phủ ~100%. Chuẩn hoá RDKit, bỏ muối —
   tên trong `drugsfda` là dạng muối ("taletrectinib adipate").
   *Không dựng pipeline UNII* — bulk file trả 403 và PubChem+ChEMBL đã đủ.
4. Lọc bỏ peptide/oligo/dược chất phóng xạ/chất cản quang/hỗn hợp theo MW và nguyên tố
   → ~270–285. *(Tỷ lệ 12–18% này là ước lượng từ quan sát tên, chưa phân loại hệ thống.)*
5. Trích xuất bằng **LLM trên mục PK**, xuất có cấu trúc như mục 3. Loại mọi
   `is_apparent = true` khỏi tập CL/VDss.
6. Quy đổi: t½ → giờ; PPB → %gắn, **giữ nguyên trạng thái kiểm duyệt**; Vd L → L/kg ở
   70 kg và **ghi rõ giả định trong bài**, lan truyền thành dải bất định chứ không giấu.
7. Khử trùng với SMILES huấn luyện TDC (InChIKey + cả mức scaffold) **và** với
   `FDA_v2.csv` / `External_test_315.csv` của PKSmart, để chữ "external" đứng vững trước
   phản biện.
8. Riêng 35–45 chất có CL/VDss: tải hồ sơ thẩm định và **kiểm tay**. accessdata **chặn
   request không có User-Agent trình duyệt** (trả `excessive-requests-apology.html`);
   dùng UA Chrome + giãn ≥3 s. Tên file theo thời kỳ: `IntegratedR`/`MultidisciplineR`
   sau ~2020, `ClinPharmR` trước đó.

## 7. Công trình liên quan cần đọc toàn văn trước khi trích dẫn

- **Koleiev et al., JCIM 2026, 66(14), 8045–8058**, doi `10.1021/acs.jcim.6c00819` —
  kiểm toán top-3 của cả 22 bảng ADMET. Theo tóm tắt tìm được: chỉ CaliciBoost, MapLight,
  MapLight+GNN qua được mọi kiểm tra; truy ra rò rỉ dữ liệu ở MiniMol, GradientBoost,
  XGBoost; khuyến nghị test ẩn và chuẩn hoá môi trường. Bản preprint bioRxiv doi
  `10.64898/2026.02.26.708193`.
  **⚠ Chưa đọc được toàn văn (bioRxiv trả 429).** Chưa xác minh kết luận cho riêng năm
  endpoint của ta, cũng chưa rõ họ có phát hành code/data không. **Phải đọc trước khi
  trích dẫn.** Lưu ý: hai bài này lập luận đúng theo hướng ta đang làm → móc nối tốt cho
  phần Introduction.
- **"What Do Oral Drugs Really Look Like?"**, J Med Chem, doi
  `10.1021/acs.jmedchem.5c02863` — 104 thuốc uống phân tử nhỏ duyệt 2020–2024, có Vd,
  CL, t½, PPB. Mẫu curation tốt nhất để tham chiếu. Nó độc lập xác nhận mọi phát hiện về
  đơn vị ở mục 3, kể cả việc trộn lẫn giá trị biểu kiến mà **không** đối chiếu lại.
- **PKSmart**, github `srijitseal/PKSmart` — tái dùng được 42 dòng cho 2021–2022 nhưng
  **phải kiểm lại từng giá trị** (mục 3).
- **PK-DataBase**, github `ClickFF/PK-DataBase` — 578 mục từ phụ lục Goodman & Gilman.
  **Thuốc cũ**, không dùng làm holdout được; chỉ dùng đối chiếu chéo.

## 8. Chưa khảo sát

EMA EPAR (không thấy API, và dính cùng vấn đề CL/F), PK-DB (pk-db.com), PharmGKB,
Open Systems Pharmacology. Không chặn việc hiện tại.

## Việc tiếp theo

- [ ] Đọc toàn văn Koleiev et al. trước khi trích dẫn.
- [ ] Dựng bước 1–3 của đường ống, báo cáo n thật thay cho ước lượng.
- [ ] Kiểm lại 42 dòng của PKSmart, ghi lại số nào bị `/F` nhiễm.
