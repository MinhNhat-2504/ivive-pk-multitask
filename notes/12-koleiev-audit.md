# Finding 12 — Koleiev et al. đọc toàn văn: ba thứ buộc phải sửa trong bài

Ngày: 2026-10-17. Đã đọc toàn văn **cả hai bản**: JCIM 2026, 66(14), 8045–8058
(doi `10.1021/acs.jcim.6c00819`, **open access**, PMC13417885) và preprint bioRxiv
(doi `10.64898/2026.02.26.708193`). Bản JCIM có nhiều nội dung preprint không có —
**trích dẫn bản JCIM**.

Tác giả đều là nhân viên Receptor.AI Inc.; ba người giữ cổ phần (đã khai báo xung đột
lợi ích).

## 1. ⚠ SỬA: con số "đầu bảng" mình đang dùng là của một mô hình **không tái lập được**

Mô hình đứng đầu `Clearance_Hepatocyte_AZ` **và** `Half_Life_Obach` là **CFA**. Koleiev
et al. loại nó ở **giai đoạn 1**: *"not accessible due to the broken links on both the
TDS website and the corresponding paper."*

Nghĩa là **0,536** và **0,576** — hai con số mình đã dùng xuyên suốt làm mốc "đầu bảng" —
đến từ một mô hình **không ai tái lập được**.

Mốc tái lập được là **MapLight+GNN**. Số họ tái lập (bảng 3 JCIM):

| Endpoint | TDC công bố | **Koleiev tái lập** | Baseline của ta |
|---|---|---|---|
| `cl_hep` | 0,498 ± 0,009 | **0,500 ± 0,012** | **0,490 ± 0,054** |
| `cl_mic` | 0,630 ± 0,010 | **0,635 ± 0,007** | 0,577 ± 0,051 |
| `vdss` | 0,713 ± 0,007 | **0,716 ± 0,006** | 0,686 ± 0,045 |
| `half_life` | 0,557 ± 0,034 | **0,563 ± 0,009** | 0,452 ± 0,087 |
| `ppbr` (MAE) | 7,526 ± 0,106 | **7,573 ± 0,163** | 8,62 ± 0,75 |

**Baseline của ta trên `cl_hep` là 0,490 so với mốc tái lập được 0,500** — cách nhau
0,010, chứ không phải 0,046 so với CFA. Mục 4.1 của bài hiện viết "mọi giá trị đầu bảng
cách trung bình của ta 3–8 lần sai số chuẩn" — **phải viết lại**.

## 2. ⚠ SỬA: phát biểu về phương sai giao thức của mình sai ở chỗ nào

Giao thức thật của `admet_group`: **test set cố định** (243 phân tử, hash không đổi), seed
chỉ chia lại train/valid (848/122). Đã tự kiểm chứng.

Mình đo bằng `get_split(method='scaffold', seed=s)` — **chia lại toàn bộ, kể cả test set**.
Nên con số 0,045–0,087 của mình **không phải phương sai của giao thức bảng xếp hạng**.
Đang đo lại bằng giao thức thật (`src/eval/leaderboard_protocol.py`).

Koleiev báo SD **0,005–0,016** trên năm endpoint này, và repo của họ ghi rõ dùng
*"the standard 5-seed protocol from `tdc.benchmark_group.admet_group`"*. Chênh 4–8 lần so
với số của mình — **khớp đúng với việc hai bên đo hai thành phần phương sai khác nhau**.

> Phát biểu đúng phải là: **độ lệch ± mà bảng công bố chỉ gồm phương sai phân hoạch
> train/valid và tính ngẫu nhiên của bộ học, vì test set cố định. Nó không nói gì về việc
> con số đó sẽ dịch bao nhiêu trên một mẫu phân tử khác.**

Đây vẫn là một phát hiện, chỉ là phát hiện **khác** với điều mình đã viết — và nó mạnh
hơn, vì nó giải thích được tại sao hai phép đo lệch nhau 4–8 lần.

*(Báo cáo tra cứu còn nêu: trong `evaluate.py` của họ, vòng lặp seed chỉ đổi
`random_state` của bộ học, train_val giữ nguyên — tức hẹp hơn nữa. Nhưng đó là pipeline
nội bộ của họ, chưa kiểm cho các script sinh ra bảng 3. **Chưa xác minh, không viết vào
bài như sự thật.**)*

## 3. ⚠ SỬA: không được nói rò rỉ của MiniMol chạm tới endpoint của mình

Rò rỉ MiniMol là **gián tiếp, ở mức kho tiền huấn luyện**: lọc theo canonical SMILES
không xử lý đúng lập thể, nên đồng phân lập thể khác của cùng một chất vẫn sót lại (ví dụ
danazol). Nhưng:

- Phân tích so với **hợp của toàn bộ phân tử test TDC**, **không tách theo endpoint**.
- **Không định lượng** bao nhiêu phân tử bị nhiễm.
- MiniMol bị loại ở giai đoạn 2 và **không bao giờ được đánh giá lại** → **không có con số
  hiệu chỉnh nào** cho bất kỳ endpoint nào.
- MiniMol đứng đầu 7 endpoint, **không endpoint nào thuộc năm của mình**.

→ Phát biểu trung thực: *một cơ chế nhiễm hợp lý đã được chứng minh ở mức kho dữ liệu và
chưa bao giờ được giải quyết theo từng endpoint* — **không phải** "số của MiniMol trên
endpoint của ta bị thổi phồng".

Endpoint duy nhất của mình có kết luận rò rỉ xác nhận là **`ppbr`**, và là của
**GradientBoost**, không phải MiniMol: lỗi gọi hàm chia lại dữ liệu với seed cứng 42
trong vòng tinh chỉnh, khiến validation set chồng lấn test set. Rank 1 → 3. **Chỉ công bố
thứ hạng, không công bố giá trị đã hiệu chỉnh.**

## 4. ✅ GIỮ: họ **không** đụng tới phát hiện kiểm duyệt của mình

Grep toàn văn cả hai bản + phụ lục cho: censor, detection limit, LLOQ, limit, truncat,
saturat, µL/min, 150.0, 3.0 → **không có kết quả nào**.

Phát hiện 27–34% nhãn nằm đúng ngưỡng phát hiện (`notes/04`) **hoàn toàn chưa ai chạm tới**.

Và bản JCIM có hẳn một lời mời:

> *"A systematic Walters-style audit of all 22 TDC ADMET end points… lies outside the
> end-user scope of the present work, but we consider it one of the most valuable
> directions for follow-up and hope the audit reported here helps motivate it."*

## 5. ✅ GIỮ: họ **không** đo phương sai giao thức, và tự nói là không làm kiểm định

Bản JCIM thêm đoạn miễn trừ:

> *"We did not apply formal statistical significance tests (e.g., Friedman/Nemenyi-style
> posthoc tests on per-fold metrics…) to the reproduced leaderboard."*

Và dẫn Kamuntavičius et al. rằng *"many of the small performance gaps separating
top-ranked TDC ADMET entries are not statistically distinguishable"*, kết luận
*"leaderboard ranking alone… should not be treated as a definitive ordering of methods."*

→ Đóng góp của mình nằm đúng vào khe họ đã nêu tên và từ chối lấp.

## 6. ✅ Bảng của chính họ ủng hộ luận điểm nhiễu của mình

Trong bảng 3, `VDss` và `Half_Life`: **cùng hai mô hình đổi chỗ cho nhau theo hai chiều
ngược nhau**, trên chênh lệch 0,001–0,009 — tức hoán vị hạng 1↔2 và 2↔3 nằm gọn trong độ
lệch chuẩn họ tự ghi. Trích dẫn được trực tiếp.

Và thí nghiệm overfit cố ý của họ: `Clearance_Microsome` từ hạng 14 lên **hạng 1**, mức
tăng **+0,077 Spearman** — nằm trong dải nhiễu mình đo. **Nhưng các số ở bảng S3/S4 là
chạy một lần, không có thanh sai số nào.** Tức luận điểm overfit trung tâm của họ không có
kiểm soát phương sai — đúng chỗ kết quả của mình bổ sung được.

## 7. Ba mô hình qua được, và "qua được" nghĩa là gì

**CaliciBoost, MapLight, MapLight+GNN** — 3/10. Nhưng 5/7 bị loại ở giai đoạn 1 là **lỗi
đóng gói phần mềm, không phải lỗi khoa học**: họ cố ý không liên hệ tác giả, không mở
issue, coi *"ability of the public artifacts alone"* là một phần của thứ đang được đánh
giá. Nếu trích dẫn con số "chỉ 3/10 qua được" thì **phải nói rõ điều này**.

CaliciBoost chỉ tồn tại cho Caco-2 → không nói gì về năm endpoint của mình.

## 8. Code và dữ liệu

Chỉ bản JCIM có: `github.com/receptor-ai/tdc-admet-bench` (có thật, cập nhật 05/2026,
không có license) và Zenodo `20180944` (**chưa xác minh được — API trả 403**).

Hai điểm đáng lưu ý: README thừa nhận *"The Mol2Vec rows from Section 3 rely on an
internal Receptor.AI package that is not public and cannot be reproduced from this repo"*;
và snapshot dữ liệu trong repo là **2026-03-24** trong khi bảng so sánh trong bài là
**11/2025** — lệch phiên bản, trong khi chính họ lấy việc trôi phiên bản dữ liệu làm luận
điểm trung tâm.

## 9. Chưa ai phản biện

JCIM 5 trích dẫn, preprint 3, đều ngẫu nhiên. Không có phản biện, không có phản hồi từ
TDC, không có issue nào trên `mims-harvard/TDC` nhắc tới. Không thấy thay đổi nào ở TDC.

Công trình liên quan đáng theo: Schuh, Daniluk & Sieber, *Chem Sci* 2026 (PMC13425739) —
kiểm toán 50+ cấu hình dataset gồm cả TDC. Và một bài **AssayInspector** (*J Cheminform*)
được nêu là báo cáo lệch giữa TDC và dữ liệu ADME chuẩn **riêng cho half-life và
clearance** — nếu đúng thì trúng endpoint của mình. **Cả hai chưa đọc.**

## Việc tiếp theo

- [ ] Chờ `src/eval/leaderboard_protocol.py` xong → sửa `notes/02` và Mục 4.
- [ ] Sửa mọi chỗ trong bài dùng 0,536 / 0,576 làm mốc "đầu bảng".
- [ ] Đọc Schuh et al. và AssayInspector.
