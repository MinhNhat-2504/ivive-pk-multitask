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

## 2. Kết quả

| % nhãn | n | đơn nhiệm | ràng buộc | **Δ** | KTC 95% | thắng | p |
|---|---|---|---|---|---|---|---|
| 25% | 79 | 0,0923 | 0,2405 | **+0,1482** | [+0,092, +0,208] | 20/25 | **<0,0001** |
| 50% | 164 | 0,2204 | 0,2938 | **+0,0734** | [+0,032, +0,114] | 18/25 | **0,0028** |
| 75% | 242 | 0,2248 | 0,3204 | **+0,0956** | [+0,064, +0,126] | 23/25 | **<0,0001** |
| 100% | 325 | 0,2739 | 0,3790 | **+0,1051** | [+0,068, +0,144] | 23/25 | **<0,0001** |

**Ràng buộc giúp ở MỌI mức nhãn, kể cả 100%** — nơi trên dữ liệu TDC gộp loài nó **có hại**
(−0,0175). Và mức giúp lớn gấp 2–3 lần.

Hợp lý: ở 100% nhãn đây vẫn chỉ là **325 phân tử huấn luyện**, trong khi tập gộp có 970.
Theo đúng câu chuyện giảm nhãn, 325 nhãn tế bào gan người **vẫn là khan hiếm**.

## 3. ⚠ Nhưng ablation đảo chiều — và đây mới là điều quan trọng

| % nhãn | identity − đơn nhiệm | p | **Deming − identity** | thắng | p |
|---|---|---|---|---|---|
| 25% | **+0,1391** | 0,0001 | +0,0092 | **4/25** | 0,080 |
| 50% | **+0,0555** | 0,032 | +0,0179 | 14/25 | 0,114 |
| 75% | **+0,0778** | 0,0002 | +0,0178 | 15/25 | 0,191 |
| 100% | **+0,0805** | 0,0028 | +0,0246 | 17/25 | **0,015** |

Trên dữ liệu người thuần, **coupling cố định (identity) cũng giúp gần bằng coupling học
được**. Lợi thế của việc học co lại còn +0,009 đến +0,025, và **chỉ có ý nghĩa ở 100% nhãn**.

**Lưu ý một mâu thuẫn phải báo cáo:** ở 25% nhãn, trung bình chênh lệch là +0,0092 với KTC
không chứa 0, **nhưng chỉ thắng 4/25**. Nghĩa là **trung vị âm** — identity thắng ở phần
lớn seed, còn trung bình bị vài seed có chênh lệch dương rất lớn kéo lên. Ở mức này phải
báo cáo **kiểm định dấu làm chính**: identity thường tốt hơn.

## 4. Diễn giải — và nó sửa lại luận điểm của bài

Trên TDC (gộp loài), "học coupling thắng cố định coupling" rất mạnh (p < 0,0001, 46/60
ngoài miền). Trên dữ liệu người thuần, nó gần như biến mất.

Giải thích hợp lý nhất, và nó hợp dược học:

> Trên tập gộp, coupling phải bắc cầu từ **microsome người** sang **hỗn hợp tế bào gan
> người–chuột**. Một coupling cố định ở quy đổi sinh lý thuần tuý sẽ **sai hệ thống**,
> nên việc học nó hấp thụ phần lệch loài và mang lại lợi ích lớn. Trên dữ liệu người
> thuần, microsome → tế bào gan **cùng loài**, nên quy đổi sinh lý trong sách giáo khoa
> **vốn đã gần đúng**, và học thêm chỉ còn ít chỗ để cải thiện.

Nghĩa là: **một phần đáng kể của "lợi ích từ việc học phép hiệu chỉnh" ở Mục 7 là hiện
vật của việc TDC gộp loài.** Phải nói thẳng điều này.

Luận điểm còn lại, và nó vững hơn:

> **Ràng buộc IVIVE — truyền thông tin từ assay microsome sang assay tế bào gan qua quy
> đổi sinh lý — giúp đáng kể và ở mọi mức nhãn trên dữ liệu người thuần (+0,07 đến +0,15,
> p ≤ 0,003).** Việc hệ số hiệu chỉnh được học hay cố định thì ít quan trọng khi hai assay
> cùng loài; nó trở nên quan trọng khi tập dữ liệu trộn nguồn không đồng nhất.

Phát biểu này **mạnh hơn** cái cũ về mặt dược học, vì nó là lời khẳng định cho chính phép
IVIVE chứ không chỉ cho một thủ thuật khớp đường thẳng.

## 5. Giới hạn

1. 408 phân tử, test TB chỉ 82. Mọi hiệu ứng đo trên cỡ nhỏ.
2. Hiệu năng tuyệt đối thấp hơn hẳn tập gộp (0,27–0,38 so với 0,44–0,49) — đúng như kỳ
   vọng với 408 so với 1 020 phân tử.
3. Chưa chạy 60 seed. Luận điểm "ràng buộc thắng đơn nhiệm" ở đây rất mạnh (p < 0,0001,
   23/25) nên ít rủi ro, nhưng ablation thì ở ranh giới và **cần thêm seed**.
4. Chưa kiểm ngoài miền trên dữ liệu người thuần — split cụm không dùng được ở cỡ này (§1).

## Việc tiếp theo

- [ ] Chạy 60 seed cho ablation người thuần — hiện ở ranh giới.
- [ ] Viết lại Mục 7 theo phát biểu ở §4.
- [ ] Xem xét chuyển toàn bộ bài sang dữ liệu người thuần làm kết quả chính, giữ tập TDC
      gộp loài làm minh hoạ cho phần phê bình benchmark.
