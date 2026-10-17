# Finding 02 — Baseline đơn nhiệm, và độ lệch giữa seed nuốt trọn mục tiêu cải tiến

Ngày đo: 2026-10-17. Script: `src/models/baseline.py`, kết quả thô ở `results/`.
Giao thức: split scaffold của TDC, seed 1–5, LightGBM 600 cây trên 28 descriptor RDKit +
Morgan count 1024 bit. Không tinh chỉnh siêu tham số, không early stopping.

## Kết quả

| Endpoint | Chỉ số | Baseline của ta | min – max qua 5 seed | Đầu bảng TDC (theo PDF) |
|---|---|---|---|---|
| `cl_hep` | Spearman ↑ | **0,502 ± 0,055** | 0,450 – 0,582 | 0,536 |
| `cl_mic` | Spearman ↑ | **0,578 ± 0,019** | 0,548 – 0,595 | 0,630 |
| `vdss` | Spearman ↑ | **0,707 ± 0,054** | 0,633 – 0,758 | 0,713 |
| `half_life` | Spearman ↑ | **0,409 ± 0,122** | 0,264 – 0,602 | 0,576 |
| `ppbr` | MAE ↓ | **8,95 ± 0,42** | 8,46 – 9,41 | 7,44 |

## Ba điều đọc ra

### 1. Mục tiêu cải tiến nhỏ hơn độ lệch giữa seed

Mục tiêu đặt ra trong PDF là `cl_hep` 0,536 → 0,60 (+0,064) và `half_life` 0,576 → 0,63
(+0,054). Nhưng chỉ riêng việc đổi seed của phép chia scaffold đã làm điểm xê dịch:

- `cl_hep`: SD = 0,055 — **gần bằng toàn bộ mức cải tiến nhắm tới**.
- `half_life`: SD = 0,122 — **gấp hơn hai lần mức cải tiến nhắm tới**. Khoảng min–max là
  0,264–0,602, rộng 0,34, tức **gấp sáu lần** mức cải tiến.

Nghĩa là: một kết quả "cải thiện +0,05" báo cáo trên 5 seed **không phân biệt được
phương pháp với may rủi của phép chia**. Đây không còn là phỏng đoán nữa — đó là số đo
trên chính pipeline của ta.

### 2. Seed đơn lẻ của baseline thường đã chạm mức đầu bảng

`half_life` seed 2 cho **0,602**, cao hơn con số 0,576 đang đứng đầu bảng. `vdss` cho
0,707 ± 0,054, tức **ngang với 0,713 của MapLight+GNN** dù đây chỉ là LightGBM chưa tinh
chỉnh. `cl_hep` seed 3 cho 0,582, cao hơn 0,536.

Điều này củng cố nhận định trong PDF rằng boosting trên descriptor đang ngang ngửa mô
hình nền đồ thị trên nhóm endpoint này — và gợi ý rằng một phần ưu thế báo cáo của các
mô hình đầu bảng có thể đến từ lựa chọn seed/giao thức chứ không từ năng lực mô hình.
Chưa kết luận được: cần lấy độ lệch chuẩn mà chính bảng TDC công bố để so.

### 3. Vẫn còn khoảng trống tinh chỉnh thật

Không được đọc hai điểm trên thành "baseline của ta đã bằng SOTA". Trung bình của ta
thấp hơn đầu bảng ở `cl_hep` (0,502 vs 0,536), `cl_mic` (0,578 vs 0,630), `half_life`
(0,409 vs 0,576) và `ppbr` (8,95 vs 7,44). Một phần khoảng cách đó là dư địa tinh chỉnh
thật. Điều đã chứng minh được là **độ rộng của nhiễu**, không phải sự ngang bằng.

## Hệ quả cho thiết kế bài

1. **Bỏ "vượt SOTA trên TDC" khỏi luận điểm chính.** Với SD 0,055–0,122, luận điểm đó
   không bảo vệ được trước reviewer, và PDF đã cảnh báo trước.
2. **Giao thức đánh giá phải đổi.** 5 seed là không đủ. Cần ≥ 20 lần chia lặp lại, báo
   cáo khoảng tin cậy, và kiểm định ghép cặp theo từng lần chia (paired test) giữa mô
   hình có ràng buộc và không ràng buộc — so sánh ghép cặp khử được phương sai của phép
   chia, vốn là nguồn nhiễu trội ở đây.
3. **Đóng góp thứ hai của bài đã có sẵn:** định lượng phương sai giao thức của bảng xếp
   hạng ADMET. Đây là kết quả tự thân, và nó biến đúng điểm yếu mà reviewer sẽ chỉ ra
   thành một phần kết quả.
4. **Luận điểm chính dời sang:** ràng buộc IVIVE mang lại *tính nhất quán cơ chế* và
   *khả năng khái quát ra miền ngoài* — đo trên tập external theo thời gian — chứ không
   phải vài phần trăm trên bảng.

## Việc tiếp theo

- [ ] Chạy lại baseline với ≥ 20 lần chia lặp lại để có khoảng tin cậy thật.
- [ ] Lấy độ lệch chuẩn công bố của bảng TDC để so trực tiếp.
- [ ] Ước lượng trần nhiễu thực nghiệm từ các lần đo lặp trong dữ liệu thô.
- [ ] Dựng external set theo thời gian từ Drugs@FDA (thuốc duyệt sau 2016).
