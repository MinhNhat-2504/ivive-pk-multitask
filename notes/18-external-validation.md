# Finding 18 — Validation ngoài thật: thứ hạng chuyển giao được, hiệu chuẩn thì không

Ngày: 2026-10-18. Script: `src/experiments/external_validation.py`.
Huấn luyện trên `train_val` của TDC, chấm trên thuốc **duyệt từ 2017**, đã kiểm tay và
rời cả benchmark lẫn PKSmart. 10 seed.

## 1. Kết quả tổng

| Endpoint | n ngoài | benchmark Spearman | **ngoài Spearman** | **chênh** |
|---|---|---|---|---|
| `half_life` | 76 (45 không kiểm duyệt) | 0,492 ± 0,011 | **0,334 ± 0,029** | **+0,158** [+0,139, +0,175] |
| `ppb` | 128 (71 không kiểm duyệt) | 0,617 ± 0,013 | **0,689 ± 0,022** | **−0,072** (tốt lên) |

Sai số:

| Endpoint | | benchmark | ngoài |
|---|---|---|---|
| `half_life` | MAE giờ thô | 7,85 | **31,80** (gấp 4) |
| `half_life` | **sai số bội trung vị** | **1,93×** | **2,31×** |
| `ppb` | MAE điểm % | 7,83 | **11,30** (tệ 44%) |
| `ppb` | sai số bội trung vị | 1,05× | 1,07× |

> **MAE giờ thô của `half_life` phóng đại mức giảm.** Nó bị đuôi dài chi phối (max 559
> giờ). Sai số bội — cách y văn dược động học báo cáo — cho **1,93× lên 2,31×**, tức giảm
> thật nhưng khiêm tốn hơn nhiều. **Phải dẫn bằng sai số bội**, không phải MAE giờ.

## 2. Hai endpoint, hai câu chuyện khác hẳn

**`half_life` khái quát kém ra ngoài thời gian.** Mất **0,158 Spearman**, KTC không chứa
0, và sai số bội từ 1,93× lên 2,31×. Trên phần không kiểm duyệt còn tệ hơn (0,275). Thuốc
duyệt sau 2017 khác đủ nhiều để mô hình huấn luyện trên Obach (biên soạn 2008) hụt rõ rệt.

**`ppb` thì ngược: thứ hạng TỐT LÊN (0,617 → 0,689) nhưng MAE TỆ ĐI 44%.**

Đó là dấu hiệu kinh điển của **lệch hiệu chuẩn**: mô hình xếp hạng các chất mới đúng,
nhưng sai về giá trị tuyệt đối. Hợp lý — PPBR_AZ là assay nội bộ AstraZeneca còn nhãn
thuốc lấy giá trị của nhà tài trợ, nên **lệch giữa hai assay** đúng như khảo sát đã dự
đoán. Thứ hạng miễn nhiễm với lệch hệ thống; MAE thì không.

> **Thứ hạng chuyển giao được, hiệu chuẩn thì không.**

## 3. ⚠ Và đây là phát hiện quan trọng nhất — nó chống lại chính phương pháp của mình

MAE nhận thức kiểm duyệt **gần như bằng y hệt** MAE ngây thơ (11,30 so với 11,31). Nghĩa
là mô hình **gần như không bao giờ rơi đúng phía của ngưỡng**. Kiểm:

| Kiểm duyệt | n | giá trị thật TB | dự đoán TB | lệch | **đúng phía** |
|---|---|---|---|---|---|
| `gt` (">99%") | 25 | 98,6 | 92,3 | **−6,28** | **4%** |
| `lt` ("<10%") | 3 | 8,7 | 49,6 | **+40,98** | **0%** |

Hiệu chuẩn theo vùng cho thấy **co về trung bình kinh điển**:

| Vùng PPB | n | thật | dự đoán | lệch |
|---|---|---|---|---|
| < 50 | 13 | 21,3 | 57,3 | **+35,9** |
| 50–90 | 34 | 72,3 | 83,0 | +10,7 |
| 90–95 | 17 | 92,2 | 84,8 | −7,4 |
| 95–99 | 27 | 96,8 | 91,8 | −4,9 |
| **≥ 99** | 37 | 99,3 | 92,9 | **−6,4** |

Mô hình nén mọi thứ về khoảng ~85%.

### Hệ quả cho ràng buộc IVIVE — định lượng được

Trong mô hình của mình, **fu = 1 − PPB/100** là biến trung gian của ràng buộc well-stirred.
Với 37 chất ngoài có PPB ≥ 99%:

| | trung vị |
|---|---|
| fu thật | **0,0100** |
| fu dự đoán | **0,0617** |
| **hệ số sai** | **8,3 lần** |

> **Trên dữ liệu thật sự ngoài miền, phần thuốc tự do — biến trung gian của chuỗi IVIVE —
> sai trung vị 8,3 lần với các chất gắn protein cao.** Và chất gắn protein cao chiếm
> 29% tập ngoài (37/128).

Vì CL_h trong mô hình well-stirred tỷ lệ gần tuyến tính với fu khi fu·CL_int ≪ Q_h, sai
số 8,3 lần ở fu truyền gần như nguyên vẹn sang thải trừ gan dự đoán.

**Đây là giới hạn nghiêm trọng của chính phương pháp mình đề xuất, và nó chỉ lộ ra khi có
tập ngoài thật.** Không có nó thì không ai biết.

## 4. Phải viết gì vào bài

1. Chuỗi IVIVE **không thể dùng dự đoán PPB làm fu** cho chất gắn cao mà không hiệu chuẩn
   lại. Ràng buộc nên dùng **nhãn PPB khi có**, và với chất không có nhãn thì phải truyền
   bất định của fu, không truyền điểm ước lượng.
2. Nén về trung bình ở hai đuôi **chính là thứ kiểm duyệt gây ra**: mô hình học rằng
   "3" và "150", "99,95%" là giá trị thường gặp và kéo mọi thứ vào giữa. Nối trực tiếp
   với `notes/04`.
3. `half_life` mất 0,158 ra ngoài thời gian — con số này nên đứng cạnh mọi điểm benchmark
   trong bài, vì nó là mức giảm thật khi gặp thuốc mới.

## 5. Giới hạn

1. 10 seed. Đủ cho ước lượng mức giảm (KTC chặt) nhưng không nhiều.
2. `ppb` n=128, `half_life` n=76 — nhỏ.
3. Chỉ baseline đơn nhiệm. **Chưa chạy mô hình có ràng buộc trên tập ngoài** vì ràng buộc
   hiện đặt trên cặp clearance, mà `notes/05` chứng minh không lấy được từ hồ sơ FDA.
4. Một phân tử `ppb` ngoài bị loại vì có mặt trong tập huấn luyện TDC (n 129 → 128) —
   kiểm tra rò rỉ hoạt động đúng.

## Việc tiếp theo

- [x] ~~Sai số bội cho `half_life`~~ → đã thêm: 1,93× → 2,31×.
- [ ] Thử hiệu chuẩn lại PPB trên tập ngoài, xem fu có cứu được không.
