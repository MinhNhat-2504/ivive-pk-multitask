# Finding 11 — Trích xuất tham số PK từ 194 nhãn, và độ chính xác đo được

Ngày: 2026-10-17. Module: `src/data/pk_extract.py`. Kết quả:
`data/external/pk_extractions.csv`.

Đây là bước 5 của đường ống `notes/05`, chạy trên 194 mục PK sạch về thời gian từ
`notes/05 §B4`.

## 1. Kết quả

| Tham số | n chất | không biểu kiến | Dùng được? |
|---|---|---|---|
| **`ppb`** (% gắn) | **164** | 164 | ✅ chính |
| **`half_life`** (giờ) | **143** | 96 | ✅ chính |
| `clearance` | 92 | **21** | ⚠ quá ít |
| `vd` | 6 | 3 | ❌ |

Trạng thái kiểm duyệt của `ppb`: 30 giá trị `>`, 4 giá trị `<`, 38 "xấp xỉ". **Phải xử lý
như dữ liệu kiểm duyệt**, không parse thành float trần.

So với dự kiến ở `notes/05 §B5` (~151 cho `ppb`, ~118 cho `half_life`): `ppb` **cao hơn**
(164), `half_life` **cao hơn** (143). Nhưng con số thô không nói lên gì nếu chưa đo độ
chính xác.

## 2. Loại bỏ theo ngữ cảnh — đúng cảnh báo của khảo sát

Khảo sát cảnh báo regex "clearance" cho 41% nhưng phần lớn là **thanh thải creatinin**.
Bộ lọc ngữ cảnh bắt được:

| Lý do loại | n |
|---|---|
| creatinine | 42 |
| renal clearance | 14 |
| cockcroft | 5 |
| dialysis | 4 |
| clcr | 4 |
| renal impairment | 4 |
| egfr | 3 |
| hepatic impairment | 2 |

Tổng 78 giá trị bị loại. Nếu không loại, `clearance` sẽ báo 92 chất thay vì **21 chất có
CL toàn thân thật** — đúng tỷ lệ 41% → 16% mà khảo sát đã đo độc lập.

**Kết luận: 21 chất là quá ít cho validation CL.** Nhánh CL/Vd của tập external coi như
không dùng được, trùng với kết luận `notes/05 §2` rằng nhánh in vitro cũng không lấy được.
Tập external **chỉ còn `ppb` và `half_life`**.

## 3. Độ chính xác — đo bằng cách đọc tay, ba vòng

Không tin bộ trích xuất của chính mình. Mỗi vòng: lấy mẫu ngẫu nhiên, đọc từng câu nguồn,
đối chiếu giá trị.

| Vòng | seed mẫu | đúng | sai | độ chính xác |
|---|---|---|---|---|
| 1 | 7 | 15/22 | 5 sai + 2 sai một phần | **68%** |
| 2 | 31 | 18/24 | 3 sai + 3 không kiểm được | **75%** |
| 3 | 99 | 21/24 | 3 | **87,5%** |

Mỗi vòng dùng **mẫu mới**, không phải mẫu đã sửa theo — nếu không thì chỉ là tự lừa.

### Các lỗi đã tìm ra và sửa

1. **Đòi chữ "binding" liền kề con số** (vòng 1). Thực tế luôn có tên thuốc chen giữa:
   *"protein binding **of neratinib was greater than** 99%"*. Đây là lỗi làm độ phủ `ppb`
   chỉ còn 52/194. Sửa: bộ lọc ngữ nghĩa đã ở mức câu rồi, nên bỏ ràng buộc liền kề.
   → 52 lên 171 chất.
2. **Câu bảng** (vòng 1, 2). Mục PK chứa bảng bị dàn phẳng thành văn bản; một "câu" bảng
   có thể nhắc "binding" ở ô này và chứa phần trăm lạc ở ô khác. Sửa: bắt buộc phần trăm
   nằm trong 120 ký tự của từ "bind/bound" (60 nếu câu dày số).
3. **Gán nhầm thuốc trong sản phẩm phối hợp** (vòng 2). `avutometinib` lấy giá trị 90%
   của `defactinib`. Sửa: gắn cờ khi câu nhắc tên thuốc khác. → giờ ra đúng 99%.
4. **Gán nhãn của chất chuyển hoá cho thuốc mẹ** (vòng 3). `sofpironium` lấy
   *"the major sofpironium **metabolite** (BBI-4010) had plasma protein binding around
   2,3–3,7%"*. Sửa: loại câu có chủ ngữ là chất chuyển hoá, nhưng **giữ** câu nhắc cả hai
   ("X là 96% và N-desmethyl X là 95%") vì giá trị đầu là của thuốc mẹ.

### Lỗi còn lại ở vòng 3 — chưa sửa

- **Câu trộn thức ăn**: *"no effect of admixing SOLOSEC with pudding and yogurt..."*
- **Lệch cột bảng**: `maribavir` lấy 24,9 — đó là **Vss tính bằng L**, còn "% bound" là ô
  kế tiếp trong bảng.

Hai kiểu này khó loại bằng quy tắc chung mà không hy sinh độ phủ.

## 4. Phát biểu trung thực về bộ dữ liệu này

> Bộ trích xuất đạt **độ chính xác ~87,5%** trên mẫu đọc tay 24 giá trị `ppb`. Nó sinh ra
> **ứng viên**, không phải bộ dữ liệu cuối. Trước khi dùng cho validation công bố, 164
> giá trị cần **kiểm tay** — ước tính vài giờ công cho một người.

Không được báo cáo kết quả external trong bài khi chưa làm bước kiểm tay đó. Tỷ lệ sai
12,5% là quá cao cho một bộ dữ liệu dùng để phản bác hay xác nhận mô hình.

Đây cũng chính là chỗ công trình liên quan đã sai: PKSmart công bố `Vz/F` chia 70 kg dán
nhãn VDss (`notes/05 §3`) — tức bỏ qua đúng bước kiểm này.

## 5. Giới hạn

1. Trích xuất dựa trên **quy tắc có nhận thức ngữ cảnh**, không phải mô hình ngôn ngữ.
   `notes/05` khuyến nghị dùng LLM; cách này rẻ hơn và đo được, nhưng trần độ chính xác
   thấp hơn.
2. Độ chính xác 87,5% đo trên `ppb`. **Chưa đo cho `half_life`** — cần một vòng đọc tay
   riêng.
3. `half_life`: 96/143 không biểu kiến. Obach là t½ tĩnh mạch còn nhãn phần lớn là đường
   uống; dưới động học flip-flop thì t½ uống phản ánh hấp thu.
4. PPBR_AZ là assay nội bộ AstraZeneca, nhãn lấy từ nhà tài trợ → **chờ lệch hệ thống
   giữa hai assay**; phải báo cáo tương quan hạng song song với MAE.

## Việc tiếp theo

- [ ] **Kiểm tay 164 giá trị `ppb`** trước khi dùng cho bất kỳ kết quả nào trong bài.
- [ ] Một vòng đọc tay riêng để đo độ chính xác `half_life`.
- [ ] Khử trùng với `FDA_v2.csv` / `External_test_315.csv` của PKSmart.
