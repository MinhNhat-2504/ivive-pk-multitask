# Finding 10 — Ràng buộc giúp nhiều hơn khi phải ngoại suy, và đảo dấu ở mức nhãn đủ

Ngày đo: 2026-10-17. Scripts: `src/data/ood_split.py`, `src/experiments/ood_generalisation.py`.
Kết quả: `results/ood_*.csv`.

Đây là luận điểm thứ ba mà `notes/06` nêu và chưa trả: ràng buộc cơ chế giữ dự đoán nằm
trong vùng khả dĩ **ở chỗ dữ liệu hết**.

## 1. Split của TDC không phải phép kiểm ngoại suy mạnh

Scaffold split tách theo **danh tính** khung Bemis–Murcko, không theo **khoảng cách**.
Hai khung khác nhau có thể cách nhau một nguyên tử. Đo bằng độ tương đồng Tanimoto tới
hàng xóm gần nhất trong tập huấn luyện:

| Split | NN similarity TB | trung vị | % test > 0,5 | % test > 0,7 |
|---|---|---|---|---|
| TDC scaffold, seed 1 | 0,450 | 0,385 | 34,2% | 13,2% |
| TDC scaffold, seed 2 | 0,502 | 0,486 | 47,3% | 22,2% |
| TDC scaffold, seed 3 | 0,474 | 0,423 | 41,6% | 16,9% |
| **Cụm Butina, seed 0** | **0,324** | 0,292 | 12,8% | **2,5%** |
| Cụm Butina, seed 1 | 0,334 | 0,295 | 10,6% | 4,5% |
| Cụm Butina, seed 2 | 0,295 | 0,275 | 4,4% | 1,6% |

Một phần ba tới gần một nửa phân tử test của TDC có hàng xóm huấn luyện tương đồng trên
0,5. Đây là **nội suy**, dù sổ sách scaffold nói khác.

Split cụm giữ lại nguyên chuỗi hoá học, NN similarity trung vị **0,295**. Và nó khó thật:
Spearman của baseline đơn nhiệm tụt từ **0,491** (in-domain) xuống **0,387** (OOD).

## 2. Kết quả — 25 seed, ghép cặp

### 100% nhãn

| So sánh | chênh | KTC 95% | thắng | p |
|---|---|---|---|---|
| **có ràng buộc − đơn nhiệm** | **+0,0382** | [+0,016, +0,059] | 18/25 | **0,0034** |
| identity − đơn nhiệm | +0,0156 | [−0,008, +0,039] | 17/25 | 0,18 |
| **Deming − identity** | **+0,0226** | [+0,006, +0,040] | 18/25 | **0,014** |

### 50% nhãn

| So sánh | chênh | KTC 95% | thắng | p |
|---|---|---|---|---|
| **có ràng buộc − đơn nhiệm** | **+0,0537** | [+0,014, +0,091] | 17/25 | **0,019** |
| identity − đơn nhiệm | +0,0201 | [−0,020, +0,059] | 15/25 | 0,31 |
| **Deming − identity** | **+0,0336** | [+0,018, +0,049] | **21/25** | **0,0002** |

## 3. Kết quả chính: đảo dấu ở mức nhãn đủ

| | in-domain (TDC scaffold) | **OOD (cụm)** |
|---|---|---|
| 50% nhãn | +0,0312 (p = 0,048) | **+0,0537 (p = 0,019)** |
| 100% nhãn | **−0,0175 (p = 0,063)** | **+0,0382 (p = 0,0034)** |

Ở mức nhãn đủ, trong miền thì ràng buộc **hơi có hại** — thiên lệch của nó lấn át khi dữ
liệu đã đủ. Ra ngoài miền thì **đảo dấu thành có lợi rõ rệt**. Đúng điều một tiên nghiệm
cơ chế phải làm: không đóng góp gì ở chỗ dữ liệu đã nói đủ, và đỡ cho mô hình ở chỗ dữ
liệu im lặng.

Đây là luận điểm thứ ba của `notes/06`, giờ đã có số đo.

## 4. Một bất thường ở seed 1 **không** giữ — ghi lại vì suýt viết nhầm

Chạy thử seed 1 cho `identity` thắng `Deming` **+0,066**, ngược hẳn mọi kết quả in-domain,
và tôi đã định diễn giải thành "quan hệ vật lý cố định khái quát tốt hơn coupling học
được". Trên 25 seed thì **Deming thắng identity ở cả hai mức nhãn** (+0,023 và +0,034,
p = 0,014 và 0,0002). Bất thường đó là nhiễu của một lần chia.

Cũng seed 1 đó cho baseline đơn nhiệm 0,532 — **cao hơn** cả in-domain. Trung bình 25
seed là 0,387. Seed 1 bốc được cụm dễ.

Đây chính là thứ `notes/02` đã cảnh báo, và lần này nó suýt làm tôi viết sai kết luận vào
bài.

## 5. `identity` không thắng đơn nhiệm có ý nghĩa ở đâu cả

Kể cả ngoài miền: +0,016 (p = 0,18) và +0,020 (p = 0,31). Nên **quy đổi sinh lý thuần tuý
không chuyển giao**. Việc học phép hiệu chỉnh là thứ tạo ra lợi ích, trong miền lẫn ngoài
miền — nhất quán với `notes/08`.

## 6. Kiểm soát rò rỉ chặt — và nó làm yếu kết quả ở mục 2–3

Bản ở mục 2 chỉ loại nhãn `cl_mic` của **chính** phân tử test. Chưa đủ cho một tuyên bố
ngoại suy: nhãn microsome của một chất **tương đồng gần** với phân tử test vẫn nói cho mô
hình biết về chuỗi hoá học đang bị giữ lại.

Bản chặt phân cụm **hợp** của cả hai tập (1 440 phân tử) rồi loại mọi phân tử microsome
**cùng cụm** với tập test. Còn lại **855 / 1 102** nhãn microsome dùng được.

### 100% nhãn

| So sánh | lỏng | **chặt** | KTC 95% (chặt) | thắng | p |
|---|---|---|---|---|---|
| có ràng buộc − đơn nhiệm | +0,0382 | **+0,0154** | [−0,0005, +0,031] | 19/25 | 0,019 |
| identity − đơn nhiệm | +0,0156 | −0,0178 | [−0,038, +0,002] | 8/25 | 0,12 |
| **Deming − identity** | +0,0226 | **+0,0333** | [+0,020, +0,047] | **21/25** | **0,0001** |

### 50% nhãn

| So sánh | lỏng | **chặt** | KTC 95% (chặt) | thắng | p |
|---|---|---|---|---|---|
| có ràng buộc − đơn nhiệm | +0,0537 | **+0,0190** | [−0,014, +0,052] | 15/25 | **0,37** |
| identity − đơn nhiệm | +0,0201 | −0,0185 | [−0,060, +0,016] | 11/25 | 0,60 |
| **Deming − identity** | +0,0336 | **+0,0375** | [+0,016, +0,061] | **22/25** | **0,0012** |

## 6b. Đọc lại cho đúng — hai luận điểm có số phận khác nhau

### Luận điểm "có ràng buộc thắng đơn nhiệm ngoài miền": **yếu đi nhiều**

Ở 50% nhãn nó **mất hẳn ý nghĩa thống kê** (+0,019, p = 0,37). Ở 100% nhãn nó ở ranh giới:
Wilcoxon p = 0,019 với 19/25 lần thắng, **nhưng khoảng tin cậy bootstrap chạm 0**
([−0,0005, +0,031]).

Hai thống kê này không đồng thuận, và **phải báo cáo cả hai**. Wilcoxon xét mẫu hình dấu
(19/25 thắng là mạnh); bootstrap xét trung bình, vốn bị vài seed âm lớn kéo xuống. Cách
phát biểu trung thực: *bằng chứng gợi ý nhưng chưa dứt khoát*.

Nghĩa là: **một phần lợi ích thấy được ở mục 2 đến từ nhãn microsome của các chất tương
đồng gần với tập test** — tức rò rỉ một phần, chính thứ bản chặt sinh ra để loại.

### Luận điểm "học coupling thắng cố định coupling": **vững, và còn mạnh lên**

+0,0333 và +0,0375, thắng 21/25 và 22/25, p = 0,0001 và 0,0012 — **mạnh hơn** bản lỏng.

Lý do mạnh hơn và đáng tin hơn: hai arm này thấy **đúng cùng một dữ liệu**, nên mọi rò rỉ
là **đồng pha và triệt tiêu**. Đây là so sánh nội bộ, miễn nhiễm với đúng vấn đề vừa làm
yếu luận điểm kia.

### Hệ quả cho bài

**Con số đưa vào bài phải lấy từ bản chặt**, và luận điểm chính ngoài miền phải chuyển từ
"ràng buộc thắng đơn nhiệm" sang **"học phép hiệu chỉnh thắng cố định nó"** — là thứ vừa
đo chắc chắn hơn, vừa đúng là điểm khác biệt với PEMAL.

Điều tôi đã báo cáo ở lượt trước (+0,0382, p = 0,0034) là **số của bản lỏng và không
được dùng**.

## 7. Giới hạn

1. Cụm Butina ở cutoff 0,65 — chưa kiểm độ nhạy theo cutoff.
2. Chỉ `cl_hep` làm mục tiêu.
3. Split cụm khó hơn nhưng không phải ngoại suy thời gian; tập FDA (`notes/05`) mới cho
   điều đó, và chỉ cho nhánh in vivo.

## Việc tiếp theo

- [ ] Độ nhạy theo cutoff phân cụm — hiện cố định 0,65.
- [ ] Tăng số seed cho luận điểm "ràng buộc thắng đơn nhiệm ngoài miền": ở ranh giới ý
      nghĩa, 25 seed chưa phân giải nổi.
