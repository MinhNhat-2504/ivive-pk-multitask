# 20 — Nén dự đoán PPB là tính chất của bài toán, và cái giá của nó trong CL_h

Hai việc chặn paper C, làm xong cả hai.

Paper A báo rằng một giá trị gắn protein dự đoán cho ra `fu` sai trung vị 8,3 lần
ở nhóm gắn trên 99%. Hai câu hỏi phản biện hiển nhiên: (1) đó là tật của
LightGBM hay của bài toán, và (2) 8,3 lần trong `fu` thì còn lại bao nhiêu trong
`CL_h`, khi mô hình well-stirred bão hoà về `Q_h`.

## 1. Sáu họ mô hình, cùng dữ liệu và cùng đặc trưng

`src/experiments/ppb_compression.py`. Cùng `PPBR_AZ` train_val, cùng bộ mô tả,
không tinh chỉnh họ nào — vì tinh chỉnh từng họ sẽ biến so sánh thành so sánh về
tinh chỉnh. Đo trên 128 chất ngoài miền đã kiểm tay (37 chất gắn ≥ 99%).

| Mô hình | Spearman | Biên độ dự đoán | Lệch ở <50% | Lệch ở >99% | Sai số `fu` |
|---|---|---|---|---|---|
| ridge (tuyến tính) | 0,461 | 46,5 | +47,4 | −7,9 | 8,2× |
| SVR (nhân RBF) | 0,554 | 24,3 | +60,5 | −7,8 | 11,2× |
| k-NN (k=5) | 0,305 | 73,1 | +47,0 | −15,5 | 17,9× |
| random forest | 0,662 | 59,2 | +40,7 | −7,0 | 10,1× |
| gradient boosting | 0,682 | 71,3 | +35,9 | −6,4 | 8,3× |
| mạng nơ-ron | 0,466 | 68,8 | +28,1 | −19,4 | 26,5× |

Biên độ quan sát là 93,9 điểm phần trăm. Mọi họ đều hẹp hơn, hẹp nhất là SVR
(24,3). Mọi họ đều quá dự đoán nhóm gắn yếu (+28 đến +61 điểm) và dưới dự đoán
nhóm gắn mạnh (−6,4 đến −19,4). Dấu của hai độ lệch giống nhau ở cả sáu, tức đây
là nén về giữa chứ không phải sai số ngẫu nhiên.

**Không có họ nào thoát.** Kết luận do đó nói về nhãn và về bài toán, không nói
về một bộ ước lượng.

Và điều đáng chú ý hơn: mô hình **xếp hạng tốt nhất** (gradient boosting,
ρ = 0,682) có sai số `fu` **nhỏ nhất** (8,3×), còn mạng nơ-ron xếp hạng kém
(ρ = 0,466) sai 26,5×. Nhưng quan hệ không đơn điệu — k-NN xếp hạng kém nhất
(0,305) lại sai `fu` ít hơn mạng nơ-ron. Nên phát biểu đúng là: **một thứ hạng
tốt trên bảng xếp hạng không bảo đảm `fu` dùng được**, chứ không phải xếp hạng
tốt thì `fu` tốt. Con số 8,3× trong paper A là trường hợp **tốt nhất** trong sáu,
không phải trường hợp đại diện.

### Một cột cần đọc cẩn thận

Cột độ dốc hồi quy (observed trên predicted) có trong CSV nhưng tôi bỏ khỏi
bảng trên. Nó không đáng tin ở đây: độ dốc bị suy giảm bởi sai số trong biến dự
báo, nên một mô hình chỉ đơn giản là *kém chính xác* vẫn có thể cho độ dốc dưới 1
dù biên độ dự đoán của nó hẹp hơn sự thật. Đúng như vậy với k-NN (0,485) và mạng
nơ-ron (0,738) — hai mô hình có phương sai dự đoán lớn. Hai độ lệch theo nhóm
không có chế độ hỏng đó. Đã sửa docstring của script, vốn phát biểu ngược chiều.
Cùng lỗi errors-in-variables đã gặp ở OLS/Deming trong notes/13.

## 2. Truyền `fu` sai qua well-stirred: phần tham số

`src/experiments/fu_propagation.py`. Sai số `CL_h` theo lần, khi `fu` sai 8,3
lần, tại các `CL_int` trong khoảng làm việc công bố của assay (µL/min/mg):

| Mức gắn | `fu` | CLint 3 | 10 | 30 | 50 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| 99,9% | 0,001 | 8,29 | 8,28 | 8,23 | 8,18 | 8,07 | **7,96** |
| 99,5% | 0,005 | 8,26 | 8,18 | 7,96 | 7,74 | 7,27 | 6,85 |
| 99,0% | 0,01 | 8,23 | 8,07 | 7,64 | 7,27 | 6,49 | 5,88 |
| 95,0% | 0,05 | 7,96 | 7,27 | 5,88 | 5,00 | 3,76 | **3,10** |

### Đây là điểm cốt lõi của paper C

Lý do người ta hay cho rằng sai số `fu` "tự triệt tiêu" là vì well-stirred bão
hoà: khi tỉ số chiết xuất cao, `CL_h → Q_h` bất kể `fu`. Điều đó đúng — nhưng
bão hoà cần `fu·CL_int` tiến tới `Q_h`, và **chính nhóm chất gắn mạnh thì điều
kiện đó không bao giờ thoả**.

Với chất gắn 99,9%, ở đỉnh khoảng assay (`CL_int` = 150 µL/min/mg), sau khi quy
đổi sinh lý `fu·CL_int` = 0,123 mL/min/kg, so với `Q_h` = 20,7 — tức **0,6% của
lưu lượng máu gan**. Hoàn toàn trong vùng tuyến tính. Sai số đi qua gần như
không suy giảm: 8,29 → 7,96.

Với chất gắn 95%, cùng chỗ đó `fu·CL_int` = 6,17, tức 30% của `Q_h`, và bão hoà
bắt đầu ăn: 7,96 → 3,10.

Nghĩa là **cơ chế duy nhất có thể cứu phép ngoại suy chỉ phục vụ đúng những chất
mà mô hình vốn đã dự đoán tạm được**. Ở nhóm gắn mạnh, sai số `fu` vừa *lớn
nhất* (mục 1) vừa *không được giảm nhẹ* (mục 2). Hai gradient chạy cùng chiều
chứ không bù nhau. Đó là luận điểm, và nó không hiển nhiên trước khi tính.

## 3. Truyền `fu` sai qua well-stirred: phần thực nghiệm

Đầu-cuối trên 128 chất ngoài miền: `CL_int` dự đoán từ assay microsome người
(nhánh sạch, có nhãn loài — xem notes/18), ghép với `fu` thật và `fu` dự đoán.

| Mô hình | Trung vị, tất cả | Trung vị, >99% | p90, >99% | Tỉ lệ > 2 lần |
|---|---|---|---|---|
| ridge | 3,27 | 9,96 | 46,7 | 0,64 |
| SVR | 3,03 | 10,02 | 56,3 | 0,69 |
| k-NN | 3,66 | 15,87 | 52,1 | 0,69 |
| random forest | 2,30 | 9,12 | 33,6 | 0,61 |
| gradient boosting | **2,23** | **7,21** | **27,1** | **0,55** |
| mạng nơ-ron | 4,21 | 20,48 | 116,3 | 0,67 |

Tiêu chí nghiệm thu thông dụng của IVIVE là trong vòng 2 lần. **55–69% số chất
vượt ngưỡng đó chỉ vì số hạng gắn protein**, trước khi tính bất kỳ sai số nào
của dự đoán độ thanh thải. Ở nhóm gắn trên 99%, trung vị 7,2–20,5 lần và p90 lên
tới 116 lần.

Các con số này là **chặn dưới**: `unbound_fraction` sàn `fu` ở 1e-3 (ứng với
99,9% gắn) để tránh chia cho 0, và giá trị quan sát cao nhất trong bộ là 99,92%.
Cái sàn đó cắt bớt sai số của đúng những chất gắn mạnh nhất.

## 4. Điều này cho phép và không cho phép kết luận gì

Cho phép: giá trị PPB dự đoán từ `PPBR_AZ` không dùng được làm `fu` trong một
phép tính well-stirred cho chất gắn mạnh, và đó không phải khuyết điểm của một
mô hình cụ thể.

Chưa cho phép, và paper C phải nói rõ:

- `CL_int` ở mục 3 là **dự đoán**, không phải đo. Sai số của nó trộn vào kết quả.
  Phần tham số (mục 2) không có vấn đề này và là chỗ dựa chính; mục 3 chỉ cho
  thấy độ lớn thực tế.
- 128 chất không phải mẫu đại diện của không gian thuốc; chúng là các thuốc mới
  được FDA phê duyệt sau 2017 có đủ dữ liệu trên nhãn.
- ~~Phần nén có nguyên nhân từ nhãn bị kiểm duyệt~~ — **đã trả lời trong
  notes/21, và câu trả lời là không.** `PPBR_AZ` chỉ có 9/2790 nhãn ở biên trên,
  ChEMBL `ppb_human` có 4/1614 và không một qualifier nào; kiểm duyệt không thể
  là nguyên nhân ở endpoint này. Nguyên nhân là **thang nhãn**: hàm mất mát L2
  trên `%` gắn đặt ~111 lần nhiều trọng số hơn vào một chênh lệch sai 1,25 lần
  trong `fu` so với một chênh lệch sai 20 lần. Khớp trên `log fu` hoặc `logit`
  hạ sai số `fu` ở nhóm gắn cao từ 8,71 xuống 2,35 lần *và* cải thiện cả MAE lẫn
  Spearman trên test set của benchmark, 25/25 seed.
- Chưa so với một mô hình PPB đã công bố (ví dụ của PK-SMART) — sáu họ ở trên đều
  do tôi huấn luyện trên cùng một bộ. Mệnh đề "không riêng mô hình của tôi" đã
  vững; mệnh đề "kể cả mô hình của người khác" thì chưa.

Hai gạch đầu dòng cuối là việc tiếp theo của paper C, theo thứ tự đó.
