# Finding 02 — Giao thức đánh giá của bảng TDC tự chế ra được mức "cải tiến" lớn hơn mục tiêu của dự án

Ngày đo: 2026-10-17. Scripts: `src/models/baseline.py`, `src/eval/protocol_variance.py`.
Kết quả: `results/baseline_lgbm_25_*.csv`, `results/protocol_variance.csv`.

Giao thức: split scaffold của TDC, **25 seed**, LightGBM 600 cây trên 28 descriptor RDKit
+ Morgan count 1024 bit, không tinh chỉnh siêu tham số.

## 1. Baseline đơn nhiệm, 25 seed

| Endpoint | Chỉ số | Baseline (25 seed) | min – max | KTC 95% | Đầu bảng TDC |
|---|---|---|---|---|---|
| `cl_hep` | Spearman ↑ | 0,490 ± 0,054 | 0,364 – 0,582 | [0,379 – 0,563] | 0,536 |
| `cl_mic` | Spearman ↑ | 0,577 ± 0,051 | 0,467 – 0,654 | [0,481 – 0,648] | 0,630 |
| `vdss` | Spearman ↑ | 0,686 ± 0,045 | 0,596 – 0,758 | [0,615 – 0,749] | 0,713 |
| `half_life` | Spearman ↑ | 0,452 ± 0,087 | 0,264 – 0,611 | [0,313 – 0,606] | 0,576 |
| `ppbr` | MAE ↓ | 8,62 ± 0,75 | 7,14 – 9,76 | [7,16 – 9,70] | 7,44 |

**Baseline của ta thật sự kém hơn đầu bảng.** Sai số chuẩn của trung bình 25 seed là
sd/5, nên khoảng cách tới đầu bảng là 3–8 lần sai số chuẩn ở cả năm endpoint. Có dư địa
tinh chỉnh và dư địa phương pháp thật. **Không được đọc phần dưới thành "baseline của ta
ngang SOTA".**

## 2. Năm seed không đủ để ước lượng chính độ lệch chuẩn của nó

Lấy mẫu lại 20 000 lần: rút ngẫu nhiên 5 trong 25 seed, xem một lần chạy 5 seed sẽ báo
cáo sd bằng bao nhiêu.

| Endpoint | sd thật (25 seed) | sd mà 5 seed báo cáo: trung vị | **khoảng 95%** |
|---|---|---|---|
| `cl_hep` | 0,054 | 0,050 | **[0,021 – 0,084]** |
| `cl_mic` | 0,051 | 0,050 | **[0,019 – 0,076]** |
| `vdss` | 0,045 | 0,044 | **[0,020 – 0,065]** |
| `half_life` | 0,087 | 0,085 | **[0,026 – 0,137]** |
| `ppbr` | 0,749 | 0,733 | **[0,272 – 1,092]** |

Một lần chạy 5 seed có thể báo sd **nhỏ hơn 2,5 lần** hoặc **lớn hơn 1,6 lần** sự thật.

Chính ta đã dính: lần chạy đầu với seed 1–5 cho `cl_mic` sd = **0,019**, trong khi sự
thật trên 25 seed là **0,051** — rơi đúng vào đáy của khoảng trên. Nếu dừng ở 5 seed,
ta đã kết luận `cl_mic` là endpoint ổn định nhất. Nó không phải.

**Hệ quả: dấu ± in trên bảng xếp hạng TDC, vốn tính từ 5 seed, không phải một ước lượng
đáng tin cho phương sai giao thức.**

## 3. Phép thử quyết định — cùng một mô hình, hai lần chạy, chênh nhau bao nhiêu?

Lấy **cùng một mô hình**, không đổi gì hết. Chạy hai lần theo đúng giao thức: mỗi lần 5
seed, hai lần không chung seed nào. Chênh lệch giữa hai trung bình báo cáo chính là mức
"cải tiến" mà riêng giao thức chế ra được. 20 000 lần rút.

| Endpoint | trung vị \|chênh\| | **phân vị 95 \|chênh\|** | max \|chênh\| |
|---|---|---|---|
| `cl_hep` | 0,024 | **0,067** | 0,122 |
| `cl_mic` | 0,022 | **0,063** | 0,124 |
| `vdss` | 0,020 | **0,055** | 0,104 |
| `half_life` | 0,037 | **0,105** | 0,200 |
| `ppbr` (MAE) | 0,33 | **0,91** | 1,79 |

Đọc thẳng:

> Trên `cl_hep`, cùng một mô hình chạy hai lần theo giao thức của bảng sẽ cho chênh lệch
> **≥ 0,067 Spearman trong 5% số lần**, mà không thay đổi gì về mô hình.

**Mục tiêu ban đầu của dự án là +0,064 trên `cl_hep` và +0,054 trên `half_life`.**
Cả hai đều **nhỏ hơn** phân vị 95 của chênh lệch giả (0,067 và 0,105). Nghĩa là đúng mức
cải tiến đó phát sinh thuần tuý do chọn seed, với tần suất trên 5%, khi **không làm gì cả**.

Và theo PDF, top-5 của `Clearance_Microsome` chênh nhau **0,005** — trong khi trung vị
chênh lệch giả của endpoint này là **0,022**, gấp hơn bốn lần toàn bộ khoảng cách top-5.
Thứ tự xếp hạng giữa các mô hình đầu bảng không mang thông tin.

## 4. Kết luận

Ba điều, độc lập nhau, đều đo được:

1. Phương sai giao thức lớn (sd 0,045–0,087 Spearman).
2. Năm seed không ước lượng nổi chính phương sai đó.
3. Giao thức tự chế ra được "cải tiến" lớn hơn mục tiêu mà dự án này đặt ra ban đầu.

Kết hợp với `notes/03` (trần nhiễu `cl_hep` ≈ 0,541): bảng xếp hạng ADMET trên nhóm
endpoint dược động học **không còn đo được sự khác biệt giữa các phương pháp nữa**.

## 5. Hệ quả cho thiết kế bài

1. **Bỏ hẳn "vượt SOTA trên TDC" khỏi luận điểm.** Không phải vì khó, mà vì vô nghĩa.
2. **Mọi so sánh phải ghép cặp theo từng lần chia** (paired, cùng seed cho cả hai mô
   hình). So sánh ghép cặp khử phương sai phép chia — vốn là nguồn nhiễu trội — và là
   cách duy nhất còn lại để phát hiện hiệu ứng cỡ 0,02.
3. **Dùng ≥ 25 seed**, báo cáo KTC, và báo cáo cả phân phối chênh lệch ghép cặp.
4. Mục 2 và 3 ở trên **chính là một đóng góp phương pháp luận của bài**, không phải phần
   phụ lục.

## Việc tiếp theo

- [ ] Hạ tầng so sánh ghép cặp theo seed (dùng cho mọi thí nghiệm về sau).
- [ ] Vẽ chung một hình: phân phối baseline 25 seed + trần nhiễu + điểm đầu bảng.
- [ ] Lấy độ lệch chuẩn mà chính bảng TDC công bố, để đối chiếu trực tiếp mục 2.
