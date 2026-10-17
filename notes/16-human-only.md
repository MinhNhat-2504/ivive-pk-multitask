# Finding 16 — Chạy lại trên dữ liệu người thuần: ràng buộc giúp nhiều hơn hẳn, nhưng **lý do thì khác**

Ngày: 2026-10-17. Script: `src/experiments/human_only.py`.
Dữ liệu: `CHEMBL3301372` (tế bào gan người, 408) làm mục tiêu, `CHEMBL3301370`
(microsome người, 1 102) làm nguồn. Overlap **187 phân tử**.
Split scaffold, 25 seed, kiểm soát rò rỉ ở mức cụm (còn 808/1 102 nhãn microsome).

## 1. Chẩn đoán trước: split cụm quá khó ở cỡ này

| Split trên 408 phân tử người | Spearman |
|---|---|
| Ngẫu nhiên | 0,374 ± 0,081 |
| **Cụm Butina** | **0,047 ± 0,161** |

Mô hình **học được** nhưng **không ngoại suy được** sang chuỗi hoá học mới từ 408 chất.
Split cụm không phân giải nổi hiệu ứng ràng buộc ở đây, nên dùng **scaffold split**, khớp
với Mục 7 gốc. Ghi lại như kết quả, không viết vòng.

## 2. Kết quả, 60 seed

### Ràng buộc so với đơn nhiệm — vững ở mọi mức nhãn

| % nhãn | n | đơn nhiệm | ràng buộc | **Δ** | KTC 95% | thắng | p |
|---|---|---|---|---|---|---|---|
| 25% | 81 | 0,0923 | 0,2313 | **+0,1390** | [+0,109, +0,171] | 51/60 | **<0,0001** |
| 50% | 161 | — | — | **+0,0994** | [+0,071, +0,129] | 48/60 | **<0,0001** |
| 75% | 243 | — | — | **+0,0957** | [+0,076, +0,116] | 53/60 | **<0,0001** |
| 100% | 325 | 0,2739 | 0,3848 | **+0,1109** | [+0,087, +0,135] | 56/60 | **<0,0001** |

Gần như không đổi so với 25 seed. **Ràng buộc giúp ở mọi mức nhãn, kể cả 100%** — nơi
trên tập TDC gộp loài nó có hại (−0,0175). Hợp lý: 325 phân tử huấn luyện vẫn là khan
hiếm.

### Học coupling so với cố định coupling

| % nhãn | **Δ** | KTC 95% | thắng | p | n cặp khớp coupling |
|---|---|---|---|---|---|
| 25% | +0,0060 | [+0,000, +0,013] | **7/60** | 0,062 | **15,9 ± 3,7** |
| 50% | **+0,0202** | [+0,008, +0,033] | 38/60 | **0,0025** | 30,8 ± 4,6 |
| 75% | **+0,0161** | [+0,004, +0,029] | 38/60 | **0,032** | 46,0 ± 4,7 |
| 100% | **+0,0271** | [+0,017, +0,037] | 46/60 | **<0,0001** | 61,0 ± 5,0 |

## 3. ⚠ Cơ chế của hàng 25% — và nó sửa diễn giải ở bản 25 seed

Bản trước mình viết rằng ở 25% nhãn "identity thường tốt hơn", vì chỉ thắng 4/25. **Sai.**

`fit_coupling` quay về identity khi số cặp dùng được **< 20**. Ở 25% nhãn số cặp là
**15,9 ± 3,7** (min 7), và **82% số lần rơi xuống dưới ngưỡng đó**.

> Nghĩa là ở 25% nhãn, trong phần lớn các lần chạy, **arm "Deming" chính là arm
> "identity"** — hai arm bằng nhau, chênh lệch đúng bằng 0, và quy ước `wins = (Δ > 0)`
> đếm hoà là không thắng. Tỷ lệ 7/60 **không phải bằng chứng identity tốt hơn**; nó là
> bằng chứng rằng phép học coupling **phần lớn chưa được thực hiện**.

Đây là hiện vật của **quy tắc dự phòng do chính mình đặt ra**, không phải phát hiện về
cách chọn coupling. Phải nói rõ trong bài.

## 4. Phát biểu đúng, và nó có giá trị thực hành

Khi có đủ cặp để khớp (từ ~30 trở lên), **học phép hiệu chỉnh thắng cố định nó ở cả ba
mức còn lại**: +0,020 (p = 0,0025), +0,016 (p = 0,032), +0,027 (p < 0,0001). Độ dốc khớp
được ổn định quanh **1,04–1,18**, đúng kỳ vọng tỷ lệ thuận.

> **Ràng buộc IVIVE giúp đáng kể trên dữ liệu người thuần ở mọi mức nhãn (+0,096 đến
> +0,139, p < 0,0001). Việc học phép hiệu chỉnh thay vì cố định nó thêm được +0,016 đến
> +0,027, nhưng chỉ khi có đủ khoảng 30 phân tử mang cả hai nhãn để khớp; dưới ngưỡng đó
> phương pháp tự thoái về hằng số quy đổi sinh lý.**

Đây là một hướng dẫn thực hành đo được, không chỉ là một con số.

### Và nó sửa lại luận điểm của Mục 7

Trên TDC gộp loài, overlap là 682 phân tử (bị thổi lên bởi phần chuột), nên coupling
**luôn khớp được**, và khoảng cách học–cố định lớn vì việc gộp loài làm hằng số cố định
**sai hệ thống**. Trên người thuần, cùng loài, quy đổi sách giáo khoa vốn đã gần đúng,
nên học chỉ thêm được một phần nhỏ.

**Một phần đáng kể của "lợi ích từ việc học phép hiệu chỉnh" ở Mục 7 là hiện vật của việc
TDC gộp loài.** Phải nói thẳng.

## 5. Giới hạn

1. 408 phân tử, test TB chỉ 82. Mọi hiệu ứng đo trên cỡ nhỏ.
2. Hiệu năng tuyệt đối thấp hơn hẳn tập gộp (0,27–0,38 so với 0,44–0,49) — đúng như kỳ
   vọng với 408 so với 1 020 phân tử.
3. Ngưỡng dự phòng 20 cặp là do mình đặt, chưa tinh chỉnh. Hàng 25% nhãn đo tác dụng của
   chính ngưỡng đó chứ không đo cách chọn coupling (§3).
4. Chưa kiểm ngoài miền trên dữ liệu người thuần — split cụm không dùng được ở cỡ này (§1).

## Việc tiếp theo

- [x] ~~Chạy 60 seed cho ablation người thuần~~ → **đã chạy**, kết quả ở §2–4.
- [ ] Quét ngưỡng dự phòng (hiện 20 cặp) để tách nó khỏi kết luận ở hàng 25% nhãn.
- [ ] Viết lại Mục 7 theo phát biểu ở §4.
- [ ] Xem xét chuyển toàn bộ bài sang dữ liệu người thuần làm kết quả chính, giữ tập TDC
      gộp loài làm minh hoạ cho phần phê bình benchmark.
