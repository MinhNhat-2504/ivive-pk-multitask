# Finding 09 — Kết quả giảm nhãn không phụ thuộc vào hai tham số chưa tinh chỉnh

Ngày đo: 2026-10-17. Script: `src/experiments/sensitivity.py`.
Kết quả: `results/sensitivity_*.csv`.

`pseudo_weight` và `n_rounds` được đặt một lần ở 0,3 và 2 rồi không xem lại. Nếu lợi ích
ở `notes/08` chỉ tồn tại quanh hai giá trị đó thì nó là hiện vật tinh chỉnh và Mục 7 phải
viết lại. Quét 8 cấu hình, ghép cặp theo seed, 25 seed, so với **cùng một** baseline đơn
nhiệm. Mặc định đặt ở **giữa** dải để nếu nó là đỉnh may mắn thì sẽ lộ ra.

## 25% nhãn

| pseudo_weight | n_rounds | chênh | KTC 95% | thắng | p |
|---|---|---|---|---|---|
| **0,1** | 2 | **+0,0548** | [+0,023, +0,090] | 18/25 | 0,0034 |
| 0,2 | 2 | +0,0535 | [+0,022, +0,090] | 18/25 | 0,0046 |
| **0,3 (mặc định)** | **2** | **+0,0529** | [+0,021, +0,088] | 18/25 | 0,0061 |
| 0,5 | 2 | +0,0487 | [+0,018, +0,082] | 17/25 | 0,0051 |
| 0,8 | 2 | +0,0499 | [+0,018, +0,085] | 18/25 | 0,0105 |
| 1,0 | 2 | +0,0503 | [+0,017, +0,089] | 18/25 | 0,0081 |
| 0,3 | 3 | +0,0472 | [+0,015, +0,084] | 17/25 | 0,0147 |
| 0,3 | 4 | +0,0460 | [+0,013, +0,082] | 17/25 | 0,0275 |

**Cả tám cấu hình đều có ý nghĩa thống kê** (p từ 0,003 đến 0,028). Biên độ chênh lệch
giữa cấu hình tốt nhất và kém nhất chỉ **0,009**, trên một dải `pseudo_weight` rộng **10 lần**.

## 50% nhãn

| pseudo_weight | n_rounds | chênh | KTC 95% | thắng | p |
|---|---|---|---|---|---|
| **0,1** | 2 | **+0,0336** | [+0,005, +0,064] | 17/25 | 0,016 |
| 0,2 | 2 | +0,0329 | [+0,005, +0,062] | 18/25 | 0,037 |
| **0,3 (mặc định)** | **2** | **+0,0312** | [+0,002, +0,062] | 16/25 | 0,048 |
| 0,5 | 2 | +0,0294 | [+0,001, +0,060] | 17/25 | 0,107 |
| 0,8 | 2 | +0,0262 | [−0,003, +0,057] | 17/25 | 0,120 |
| 1,0 | 2 | +0,0251 | [−0,004, +0,056] | 17/25 | 0,149 |
| 0,3 | 3 | +0,0300 | [−0,000, +0,062] | 17/25 | 0,052 |
| 0,3 | 4 | +0,0261 | [−0,005, +0,058] | 16/25 | 0,127 |

Mọi ước lượng điểm đều dương (+0,025 đến +0,034), thắng 16–18/25 ở mọi cấu hình. Một số
cấu hình mất ý nghĩa thống kê (p tới 0,149) nhưng **KTC chồng lấn gần như hoàn toàn** —
đó là vấn đề lực thống kê, không phải bất ổn định.

## Ba điều đọc ra

### 1. Không phải hiện vật tinh chỉnh

`pseudo_weight` thay đổi 10 lần (0,1 → 1,0) làm chênh lệch dịch chưa tới 0,009 ở 25%
nhãn. Kết luận của `notes/08` đứng vững.

### 2. Mặc định không phải tối ưu — và nó làm kết quả **nhẹ đi**, không nặng thêm

Ở cả hai mức nhãn, `pseudo_weight = 0,1` tốt hơn mặc định 0,3 (+0,0548 so với +0,0529;
+0,0336 so với +0,0312). Nghĩa là con số đang báo cáo **thấp hơn** mức tốt nhất tìm được.

**Nhưng vẫn giữ mặc định 0,3 làm con số chính.** Đổi sang 0,1 vì nó cho kết quả đẹp hơn
là chọn theo kết quả — đúng thứ `notes/02` đã chứng minh là nguy hiểm ở đây. Mặc định
được chọn trước khi biết kết quả, nên nó là con số trung thực để báo cáo; cả bảng quét
được công bố kèm.

### 3. Nhiều vòng lặp hơi có hại, đơn điệu

n_rounds 2 → 3 → 4 cho +0,0529 → +0,0472 → +0,0460 ở 25% nhãn, và +0,0312 → +0,0300 →
+0,0261 ở 50%. Giảm đều ở cả hai mức. Hợp lý: mỗi vòng thêm làm pseudo-target trôi xa
nhãn thật hơn — thiên lệch tự khẳng định tích luỹ. Hai vòng là đủ, và đó cũng là mặc định
đã chọn.

## Giới hạn

- Chỉ quét hai tham số. `LABEL_WEIGHT`/`PREDICTION_WEIGHT` (0,05 so với 1,0) trong phép
  chiếu chưa quét — đó là tham số thứ ba chưa tinh chỉnh.
- Chỉ quét ở 25% và 50% nhãn, là nơi có hiệu ứng. Chưa kiểm ở 100% nhãn, nơi hiệu ứng âm.
- Siêu tham số LightGBM giữ cố định suốt, chung cho mọi arm.

## Việc tiếp theo

- [ ] Quét `LABEL_WEIGHT`/`PREDICTION_WEIGHT` để khép nốt tham số thứ ba.
- [ ] Kiểm ngoài miền — luận điểm thứ ba của `notes/06`, vẫn còn nợ.
