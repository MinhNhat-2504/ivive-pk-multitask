# Finding 08 — Ràng buộc IVIVE giảm nhu cầu nhãn, và lợi ích đến từ việc *học* phép hiệu chỉnh

Ngày đo: 2026-10-17. Scripts: `src/models/constrained.py`, `src/experiments/label_efficiency.py`.
Kết quả: `results/label_efficiency_*.csv`.

Ghép cặp theo seed, 25 seed, 4 mức nhãn. Ba arm, cùng một test fold `cl_hep`:

| Arm | Thấy gì |
|---|---|
| `single_task` | chỉ phần nhãn `cl_hep` được giữ lại |
| `constrained` | thêm nhãn `cl_mic` + coupling **học bằng Deming** |
| `constrained_identity` | thêm nhãn `cl_mic` + coupling **cố định** ở slope 1, intercept 0 |

Arm thứ ba là ablation bắt buộc: nó giả định "quy đổi sinh lý là đủ", đúng kiểu PEMAL cố
định hệ số trước. Nó thấy **đúng cùng lượng dữ liệu** với arm Deming, chỉ khác ở chỗ
không học gì từ dữ liệu.

Kiểm soát rò rỉ: **139 nhãn `cl_mic` bị loại mỗi seed** vì phân tử nằm trong test fold.

## 1. Có ràng buộc so với đơn nhiệm

| % nhãn | n | đơn nhiệm | có ràng buộc | **chênh** | KTC 95% | thắng | p |
|---|---|---|---|---|---|---|---|
| 25% | 241 | 0,3247 | 0,3776 | **+0,0529** | [+0,021, +0,087] | 18/25 | **0,0061** |
| 50% | 482 | 0,3972 | 0,4284 | **+0,0312** | [+0,002, +0,062] | 16/25 | **0,048** |
| 75% | 725 | 0,4505 | 0,4560 | +0,0055 | [−0,012, +0,022] | 15/25 | 0,58 |
| 100% | 970 | 0,4911 | 0,4736 | −0,0175 | [−0,035, +0,000] | 8/25 | 0,063 |

Đơn điệu giảm, cắt 0 quanh 75%, hơi có hại ở 100%. Đúng hành vi của một tiên nghiệm có
thiên lệch nhưng mang thông tin.

**Quy ra số nhãn tiết kiệm được** (nội suy trên đường đơn nhiệm):

| % nhãn | n thật | tương đương n đơn nhiệm | **hệ số** |
|---|---|---|---|
| 25% | 241 | 417 | **1,73×** |
| 50% | 482 | 624 | **1,29×** |
| 75% | 725 | 758 | 1,05× |

## 2. Ablation — lợi ích đến từ đâu?

### 2a. Chỉ thêm dữ liệu, không học coupling (`identity` so với đơn nhiệm)

| % nhãn | chênh | KTC 95% | p |
|---|---|---|---|
| 25% | +0,0376 | [+0,005, +0,077] | 0,090 |
| 50% | +0,0079 | [−0,021, +0,039] | 0,751 |
| 75% | −0,0080 | [−0,028, +0,013] | 0,634 |
| 100% | **−0,0291** | [−0,046, −0,011] | **0,0038** |

**Thêm dữ liệu mà không học phép hiệu chỉnh không giúp có ý nghĩa ở bất kỳ mức nào, và
có hại rõ rệt khi nhãn đủ.**

### 2b. Học coupling so với cố định coupling (`Deming` so với `identity`)

Đây là so sánh **sạch nhất** — hai arm thấy đúng cùng dữ liệu, chỉ khác cách khớp coupling.

| % nhãn | chênh | KTC 95% | thắng | p |
|---|---|---|---|---|
| 25% | +0,0153 | [+0,005, +0,026] | 16/25 | **0,016** |
| 50% | **+0,0234** | [+0,014, +0,033] | **21/25** | **0,0001** |
| 75% | +0,0134 | [+0,001, +0,025] | 15/25 | 0,052 |
| 100% | +0,0116 | [+0,002, +0,022] | 14/25 | 0,114 |

**Học phép hiệu chỉnh thắng cố định nó ở cả bốn mức nhãn**, KTC chặt (±0,01) vì đây là
so sánh ghép cặp khớp nhất có thể. Mạnh nhất ở 50%: p = 0,0001, thắng 21/25.

## 3. Phân rã trung thực

Ở mức nhãn thưa (25%), **cả hai thành phần cùng đóng góp**: thêm dữ liệu cho +0,038
(chưa đủ ý nghĩa, p = 0,09) và học coupling cộng thêm +0,015 (p = 0,016). Không được nói
toàn bộ +0,053 là nhờ ràng buộc.

Ở mức nhãn đủ (100%), bức tranh đảo: thêm dữ liệu **có hại** (−0,029, p = 0,004), và học
coupling chỉ kéo lại được +0,012 nên tổng vẫn âm.

Phát biểu đúng:

> Nhãn microsome thay thế được cho nhãn hepatocyte ở tỷ lệ đo được (~1,7× khi nhãn thưa),
> **nhưng chỉ khi phép hiệu chỉnh giữa hai hệ thống được học từ dữ liệu.** Chuyển giao
> ngây thơ — giả định quy đổi sinh lý là đủ — không giúp ở đâu cả và có hại khi nhãn đủ.

Đây đúng là điểm khác PEMAL, và giờ nó có số đo chứ không phải lập luận.

## 4. Vì sao Deming quan trọng đến thế — ghi lại vì suýt bỏ sót

Bản đầu khớp coupling bằng OLS, cho hệ số dốc **0,40–0,49**. Vô lý: hai assay đo cùng một
đại lượng thì sau quy đổi phải tỷ lệ thuận. Nguyên nhân là **suy giảm hồi quy** — OLS giả
định biến x không nhiễu, nhưng `log CL_hep` có ICC ≈ 0,57 (`notes/04`), nên hệ số dốc bị
kéo về 0 theo đúng hệ số tin cậy đó.

Deming (sai số ở cả hai biến, tỷ lệ phương sai = 1) cho **0,91–1,11** — đúng kỳ vọng vật
lý. Và kết quả ở 100% nhãn đảo từ **−0,018 sang +0,005** trên seed thử.

OLS ở đây không phải lựa chọn tồi, nó là **estimator sai**.

## 5. Giới hạn — phải viết vào bài

1. Chỉ kiểm trên `cl_hep` làm mục tiêu. Chưa thử chiều ngược (`cl_mic` làm mục tiêu).
2. Việc rút nhãn là ngẫu nhiên đều, không mô phỏng cách nhãn thật sự thưa trong thực tế
   (thường thiên lệch theo loạt chất, theo dự án).
3. Arm có ràng buộc dùng ~1 011 nhãn `cl_mic`. Hệ số 1,73× là **tỷ lệ thay thế**, không
   phải "được không mất gì".
4. Chưa kiểm ngoài miền — đó là luận điểm thứ ba mà `notes/06` nêu, vẫn còn nợ.
5. `pseudo_weight = 0,3` và `n_rounds = 2` chọn chưa qua tinh chỉnh. Cần kiểm độ nhạy.

## Việc tiếp theo

- [ ] Độ nhạy theo `pseudo_weight` và `n_rounds`.
- [ ] Chiều ngược: `cl_mic` làm mục tiêu, `cl_hep` làm nguồn.
- [ ] Kiểm ngoài miền trên tập external `ppbr`/`half_life` (`notes/05`).
