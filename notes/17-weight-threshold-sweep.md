# Finding 17 — Hai tham số cuối: không cái nào chi phối kết quả

Ngày: 2026-10-18. Script: `src/experiments/weight_sweep.py`.
Kết quả: `results/weight_sweep_*.csv`. Chạy trên dữ liệu người thuần, 25 seed,
ghép cặp với cùng baseline đơn nhiệm.

`notes/09` đã quét `pseudo_weight` và số vòng lặp. Còn hai, và cả hai đều đang
**hard-code** nên trước hết phải đưa ra thành tham số.

## 1. Trọng số phép chiếu — dải 100 lần, không đổi kết luận

Phép chiếu kéo hai dự đoán clearance của một phân tử lên đường ràng buộc, có trọng số sao
cho head **đang giữ nhãn** gần như không dịch còn head **chỉ dự đoán** hấp thụ phần hiệu
chỉnh. Tỷ lệ cũ cố định ở 0,05 so với 1,0, và biện minh dựa vào độ tản trong phân tử —
thứ `notes/15` đã **rút lại**. Nên nó cần được quét, không phải được trích dẫn.

| `label_weight` | 25% nhãn | p | 100% nhãn | p |
|---|---|---|---|---|
| 0,01 | +0,1415 | 0,0001 | +0,0956 | <0,0001 |
| **0,05 (mặc định)** | **+0,1482** | <0,0001 | **+0,1051** | <0,0001 |
| 0,10 | +0,1444 | <0,0001 | +0,1081 | <0,0001 |
| 0,25 | **+0,1578** | <0,0001 | +0,1074 | <0,0001 |
| 0,50 | +0,1338 | <0,0001 | +0,0971 | <0,0001 |
| **1,00** | +0,1243 | <0,0001 | +0,0992 | <0,0001 |

Dải **100 lần**, biên độ chỉ 0,034 (25% nhãn) và 0,013 (100% nhãn), **mọi cấu hình đều có
ý nghĩa thống kê**, thắng 20–24/25.

**Đáng chú ý: `label_weight = 1,0` — tức bỏ hẳn phân biệt giữa head có nhãn và head chỉ
dự đoán — vẫn cho +0,124 và +0,099.** Nghĩa là **cơ chế trọng số không phải thứ tạo ra
kết quả**; bản thân ràng buộc mới là. Phải nói thẳng điều này: nó làm nhẹ một phần thiết
kế mà mình từng mô tả như có chủ đích.

Mặc định 0,05 không phải tối ưu (0,25 tốt hơn ở 25% nhãn), nhưng **vẫn giữ mặc định** —
chọn theo kết quả là đúng thứ `notes/02` chứng minh benchmark này trừng phạt.

## 2. Ngưỡng dự phòng — không chi phối, và xác nhận `notes/16` từ hướng khác

| `min_pairs` | % lần coupling được khớp | 25% nhãn | 100% nhãn |
|---|---|---|---|
| 5 | 100% | +0,1487 | +0,1051 |
| 10 | 96% | +0,1498 | +0,1051 |
| **20 (mặc định)** | **20%** | +0,1482 | +0,1051 |
| 30 | 0% | +0,1391 | +0,1051 |
| 50 | 0% | +0,1391 | +0,1051 |

**Ở 100% nhãn ngưỡng không có tác dụng gì** — coupling luôn khớp (61 cặp > 50), bốn chữ
số giống hệt nhau. Đúng như dự đoán.

**Ở 25% nhãn, ngưỡng điều khiển hoàn toàn việc có khớp coupling hay không** (100% → 0%),
nhưng chênh lệch kết quả chỉ **+0,1498 so với +0,1391 = 0,011**.

> Tức là ở mức nhãn khan hiếm, việc **khớp coupling đóng góp ~0,011**, còn **toàn bộ phần
> còn lại (+0,139) đến từ chính cấu trúc ràng buộc** — truyền thông tin microsome qua quy
> đổi sinh lý.

Đây là xác nhận độc lập cho `notes/16 §3`, đo theo hướng khác: ở đó ablation Deming−identity
cho +0,006 tại 25% nhãn; ở đây quét ngưỡng cho +0,011. Cùng kết luận, cùng độ lớn.

Lưu ý: mặc định 20 rơi đúng **vùng chuyển tiếp** (20% số lần khớp) — chỗ khó diễn giải
nhất. Kết quả chính ở `notes/16` dùng giá trị này, nên hàng 25% nhãn ở đó là **hỗn hợp**
giữa có khớp và không khớp. Đã nêu trong bài.

## 3. Tổng kết cả ba lần quét

| Tham số | Dải quét | Biên độ ảnh hưởng | Có chi phối không? |
|---|---|---|---|
| `pseudo_weight` (`notes/09`) | 0,1 – 1,0 (10×) | 0,009 | ❌ |
| `n_rounds` (`notes/09`) | 2 – 4 | 0,007, giảm đơn điệu | ❌ |
| `label_weight` | 0,01 – 1,0 (100×) | 0,034 | ❌ |
| `min_coupling_pairs` | 5 – 50 | 0,011 | ❌ |

**Không tham số nào trong bốn cái chi phối kết quả.** Kết luận của bài không phải hiện vật
tinh chỉnh — và điều đó giờ **đo được**, không phải khẳng định.

## 4. Quét liên hợp — loại trừ tương tác

Quét từng tham số một không thấy được tương tác: một tham số có thể vô hại khi đứng riêng
mà vẫn quan trọng khi kết hợp. Lưới đầy đủ bốn chiều là 540 cấu hình, nên **lấy mẫu ngẫu
nhiên 24 cấu hình** từ không gian đó (gồm cả mặc định), 25 seed, cùng baseline.
Script: `src/experiments/joint_sweep.py`.

| | 25% nhãn | 100% nhãn |
|---|---|---|
| delta nhỏ nhất | **+0,1263** | **+0,0667** |
| trung vị | +0,1427 | +0,0983 |
| lớn nhất | +0,1574 | +0,1129 |
| biên độ | 0,031 | 0,046 |
| **cấu hình có delta ≤ 0** | **0/24** | **0/24** |
| **cấu hình thắng < 13/25 seed** | **0/24** | **0/24** |

> **Không một cấu hình nào trong 24 mẫu cho hiệu ứng âm hoặc thua đa số.** Cấu hình kém
> nhất ở 25% nhãn vẫn cho +0,126, ở 100% nhãn vẫn cho +0,067.

Biên độ 0,031 và 0,046 **cùng cỡ** với quét từng tham số (0,034), nên **không có tương
tác nào nổi lên**. Kết luận của quét một chiều giữ nguyên khi quét liên hợp.

Một mẫu hình nhỏ: ba cấu hình kém nhất ở 100% nhãn đều có `pseudo_weight` cao (1,0, 1,0,
0,8), khớp với `notes/09` rằng trọng số pseudo-target và số vòng lặp cao thì hại.

## 5. Giới hạn

1. Quét trên 25 seed, hai mức nhãn (25% và 100%). Chưa quét ở 50/75%.
2. Quét liên hợp lấy mẫu 24/540 cấu hình. Phủ tương tác tốt hơn lưới ở cùng ngân sách,
   nhưng không vét cạn.
3. Siêu tham số LightGBM giữ cố định suốt, chung cho mọi arm.
