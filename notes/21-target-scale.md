# 21 — Nguyên nhân của nén là thang nhãn, và sửa nó không mất gì

## Vì sao Tobit là thí nghiệm sai

notes/20 kết luận việc tiếp theo là huấn luyện trên nhãn giữ qualifier để xem
kiểm duyệt có gây ra nén. Tôi đã định làm, rồi đếm nhãn ở biên trước:

- TDC `PPBR_AZ`: 9/2790 ở đúng 99,95; **không có** giá trị nào ở biên dưới 10,0
  (nhỏ nhất là 10,09).
- ChEMBL `ppb_human`: 4/1614 ở 99,95; `standard_relation` **toàn bộ** là NaN.

Bốn trên 1614 không thể gây ra nén. **Kiểm duyệt bị loại khỏi danh sách nghi
phạm cho PPB** — khác hẳn độ thanh thải, nơi 27–34% nhãn là cờ ngưỡng phát hiện
(notes/04, paper A mục 4.2). Nếu cứ chạy Tobit theo kế hoạch, tôi sẽ thu được
"không có tác dụng" và kết luận sai rằng kiểm duyệt không quan trọng — trong khi
sự thật là ở endpoint này nó chưa bao giờ hiện diện.

Việc này cũng phơi ra một overclaim trong paper A: tiêu đề Hình 3 là "Training on
censored labels compresses predictions". Phần thân bài không hề nói vậy — dòng
475 phát biểu đúng rằng hệ số khuếch đại `1/(1 − PPB/100)` là **số học**, không
phải tính chất của mô hình. Chỉ tiêu đề hình gán nhân quả. Đã sửa thành
"Predictions compress toward the interior…" và thêm một câu giới hạn phạm vi
ngay nơi giới thiệu phép nén.

## Nghi phạm còn lại: chính thang nhãn

Nhãn của benchmark là `%` gắn và mọi hàm mất mát tính trên nó (`LGBM_PARAMS`
không đặt `objective`, nên là L2 trên `%`). Đại lượng mà IVIVE tiêu thụ là
`fu = 1 − %/100`, và ánh xạ giữa hai thang kỳ dị khi `%` tiến tới 100.

Trên tập huấn luyện: `fu` trải **1798 lần**. Trong đó

- 81,9% nhãn (dưới 99% gắn) chứa 89 lần,
- **18,1% nhãn (trên 99% gắn) chứa 20 lần, gói trong đúng một điểm phần trăm.**

Số học của gradient, ai cũng kiểm được:

| Chênh lệch | Sai số trong `fu` | Bình phương sai số trên `%` |
|---|---|---|
| 99,0 so với 99,95 | **20 lần** | 0,90 |
| 50 so với 60 | 1,25 lần | 100 |

Hàm mất mát đặt **~111 lần nhiều trọng số hơn** vào một chênh lệch tương ứng
sai 1,25 lần trong `fu`, so với một chênh lệch tương ứng sai 20 lần. Mô hình
không nén vì thiếu tín hiệu; nó nén vì được yêu cầu nén.

## Thí nghiệm: ba thang mục tiêu

`src/experiments/target_scale.py`. Cùng dữ liệu, cùng đặc trưng, cùng
`LGBM_PARAMS`, chỉ đổi thang mục tiêu. 25 seed, ghép cặp theo seed, Wilcoxon.
Sàn `fu` ở 5e-4 — đúng bằng biên 99,95% nên không cắt nhãn nào.

### Trên 128 chất ngoài miền (37 chất gắn ≥99%)

| Chỉ số | `percent` | `log_fu` | `logit` |
|---|---|---|---|
| MAE (điểm %) | 11,36 | 10,58 | **10,15** |
| Spearman | 0,686 | **0,737** | 0,733 |
| Sai số `fu`, trung vị, >99% | 8,71 | **2,35** | 2,46 |
| Sai số `fu`, p90, >99% | 30,3 | **14,8** | 16,5 |
| Tỉ lệ trong 2 lần, >99% | 0,077 | **0,408** | 0,374 |
| Biên độ dự đoán (sự thật: 93,9) | 68,1 | **91,6** | 77,8 |

### Trên test set cố định của chính benchmark (n=559, 105 chất gắn ≥99%)

| Chỉ số | `percent` | `log_fu` | `logit` | `logit` so với `percent`, ghép cặp |
|---|---|---|---|---|
| MAE (chỉ số bảng xếp hạng) | 7,812 | 7,599 | **7,433** | −0,380 [−0,556, −0,175] p=6e-8, **25/25** |
| Spearman | 0,621 | 0,679 | **0,681** | +0,060 [+0,028, +0,085] p=6e-8, **25/25** |
| Sai số `fu`, trung vị, >99% | 10,56 | **3,91** | 3,95 | −6,61 [−8,50, −5,12] p=6e-8, **25/25** |

Baseline `percent` cho Spearman 0,621, khớp với 0,617 ± 0,013 mà paper A đã báo
— nên đây không phải bù nhìn tự dựng. MAE 7,81 nằm trong vùng của các mô hình
đầu bảng TDC cho `PPBR_AZ`.

## Dự đoán của tôi sai

Tôi dự đoán sẽ có **đánh đổi**: `log_fu` tốt hơn về `fu`, tệ hơn về chỉ số
benchmark, và luận điểm của paper C sẽ là "chỉ số của benchmark và yêu cầu của
IVIVE bất đồng". Ghi lại vì dự đoán đó sai hoàn toàn. **Không có đánh đổi.**
Đổi thang tốt hơn ở mọi chỉ số tôi đo, trên cả hai bộ kiểm, 25/25 seed. Luận
điểm mạnh hơn nhưng khác: không phải xung đột giữa hai mục tiêu, mà là **một
lựa chọn mặc định sai mà không ai được lợi gì từ nó**.

(Đây là lần thứ hai trong dự án tôi đoán có đánh đổi và không có — lần đầu là
Tobit ở notes/09, nhưng ở đó dự đoán sai theo chiều ngược: tôi tưởng có lợi và
không có.)

## Confound đã loại

Phép nghịch đảo của `log_fu` chặn dự đoán ở 99,95%, còn `percent` được tự do tới
100 — có thể nghi rằng lợi thế đến từ cái chặn đó. Không:

- Chỉ **6/640** dự đoán thang `percent` vượt 99,95.
- Chặn `percent` ở 99,95 thay vì 100 làm trung vị sai số `fu` đổi **0,0000**
  (8,7096 cả hai).

Và nén nhìn thấy trực tiếp: trên bộ ngoài miền, mô hình `percent` chưa bao giờ
dự đoán dưới **27,75**, dù sự thật xuống tới 6,0. `log_fu` xuống tới 2,95.

## Việc này đổi khung của paper C

Khung cũ: "giá trị PPB dự đoán không dùng được làm `fu` cho chất gắn cao" — một
lời phàn nàn, và là lời khuyên tiêu cực.

Khung mới: **hãy khớp mô hình trên thang mà phép tính hạ nguồn tiêu thụ.** Có
chẩn đoán (số học gradient), có bằng chứng nén không phụ thuộc bộ ước lượng (sáu
họ, notes/20), có định lượng hậu quả PBPK (truyền qua well-stirred, notes/20), và
có cách sửa một dòng, cải thiện cả chỉ số của bảng xếp hạng. Đó là một bài báo
cho người làm IVIVE/PBPK *và* cho người làm benchmark.

## Cách sửa đi được ở cả sáu họ

`src/experiments/scale_across_families.py`. Sáu họ × ba thang, 5 seed, mục tiêu
được chuẩn hoá trong mọi ô để dung sai của SVR và của MLP có cùng nghĩa trên cả
ba thang.

Trung vị sai số `fu`, nhóm gắn >99%:

| Họ | `percent` | `log_fu` | `logit` | Lần cải thiện |
|---|---|---|---|---|
| ridge (tuyến tính) | 20,00 | **2,96** | 3,19 | 6,8× |
| SVR (nhân RBF) | 15,32 | **4,34** | 4,91 | 3,5× |
| k-NN (k=5) | 17,90 | **9,15** | 10,26 | 2,0× |
| random forest | 10,50 | 2,57 | **2,51** | 4,2× |
| gradient boosting | 8,28 | **2,35** | 2,58 | 3,5× |
| mạng nơ-ron | 17,38 | 4,94 | **4,28** | 4,1× |

Thang `percent` cho 8,3–20,0 lần; thang tốt nhất cho 2,4–9,2 lần. **Cả sáu họ
cải thiện, 5/5 seed mỗi họ.**

Và trên test set của benchmark, cả sáu họ cải thiện **cả hai** chỉ số:

| Họ | MAE: `percent` → tốt nhất | Spearman: `percent` → tốt nhất |
|---|---|---|
| ridge | 10,07 → 8,40 | 0,423 → 0,538 |
| SVR | 9,32 → 8,15 | 0,463 → 0,578 |
| k-NN | 12,03 → 10,01 | 0,372 → 0,392 |
| random forest | 8,55 → 7,75 | 0,607 → 0,630 |
| gradient boosting | 7,81 → 7,44 | 0,617 → 0,676 |
| mạng nơ-ron | 10,11 → 8,55 | 0,410 → 0,495 |

Sáu trên sáu, ở cả ba chỉ số. Đây là mức sạch nhất mà một kết quả có thể đạt, và
nó khép mệnh đề nhân quả: triệu chứng có ở mọi họ, cách sửa cũng hiệu nghiệm ở
mọi họ, nên nguyên nhân nằm ở thứ dùng chung — thang nhãn.

### Một lỗi tôi tự gây ra khi loại confound

`SVR(epsilon=2.0)` được chọn cho thang `%` thô, mà SD của nhãn là 18,78 — tức
ống rộng 0,106 SD. Tôi chuẩn hoá mục tiêu để dung sai có cùng nghĩa trên ba
thang, nhưng giữ nguyên `epsilon=2.0`, nên ống thành 2 SD và SVR trả về gần như
hằng số trên **cả ba** thang: biên độ dự đoán 4,5 / 0,5 / 2,9 và MAE 28,9 trên
`percent`. Dấu hiệu rất rõ — MAE 28,9 khi baseline là 9,3 — nhưng nó chỉ lộ ra
khi tôi so hàng SVR với `ppb_compression.py` (Spearman 0,554 ở đó, 0,274 ở đây).
Đã sửa bằng cách suy `epsilon` từ độ lệch chuẩn của nhãn thay vì gán tay. Bài
học lặp lại: phép "chuẩn hoá cho công bằng" tự nó là một thay đổi cần kiểm.

## Còn thiếu, theo thứ tự ưu tiên
1. **Một endpoint thứ hai.** Bài học tổng quát cần ví dụ thứ hai. Độ thanh thải
   đã dùng thang log sẵn (`ENDPOINT_SPEC`), nên ứng viên tự nhiên là `VDss` hoặc
   bán thải — nơi đại lượng tiêu thụ cũng khác thang nhãn.
2. **So với một mô hình PPB đã công bố.** Vẫn còn nợ từ notes/20; sáu họ đều do
   tôi huấn luyện.
3. Kiểm xem `logit` hay `log_fu` nên là khuyến nghị. `logit` tốt hơn về MAE,
   `log_fu` tốt hơn về `fu` ở nhóm gắn cao. Chênh lệch nhỏ so với khoảng cách
   tới `percent`; có thể khuyến nghị `logit` vì nó đối xứng ở hai đuôi và không
   cần sàn.
