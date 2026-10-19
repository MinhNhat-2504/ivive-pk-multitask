# 22 — Luận điểm chính của C đã có từ 2013; và `PPBR_AZ` gộp năm loài

Phiên này đổi hai thứ lớn: khung của paper C, và một mệnh đề sai trong paper A.

## 1. Thang biến đổi không phải phát hiện của tôi

Việc đầu tiên trong danh sách là đọc Watanabe 2018 vì nhan đề nghe quá giống
vấn đề của C. Kết quả xấu hơn tôi lo.

**Zhu, Sedykh, Zhu, Liu, Tropsha, *Pharm. Res.* 2013, 30, 1790–1798.** Trích
nguyên văn phần kết luận của abstract:

> "Models built with lnKa were significantly more accurate (MAE of 6.2–10.7 %)
> than those built with %PPB (MAE of 11.9–17.6 %) **for highly bound compounds**
> both for the training and the external sets. The pseudo binding constant
> (lnKa) is more appropriate for characterizing PPB binding than conventional
> %PPB."

`lnKa = ln(phần gắn / phần không gắn)` chính là thang `logit` của tôi. Và bài
này **đã nằm trong reference list của tôi** dưới khoá `votano2` — tôi trích nó
mà chưa đọc nó. Đây là lỗi tệ nhất trong dự án tính theo thời gian mất đi.

Hai bài nữa cùng hướng:

- **Watanabe et al. 2018**, *Mol. Pharm.* 15, 5302–5311: biến đổi logarit để xử
  lý phân bố lệch, "leading to improved accuracy at lower values".
- **Han, Xia, Xia, Tetko, Wu**, *Eur. J. Pharm. Sci.* 2025 (bioRxiv 2024.07.12):
  dùng "LogIt", sai số nhỏ hơn đáng kể so với thang phần trăm cho chất gắn trên
  99%. SOTA hiện tại, công khai trên OCHEM.

Nên mục 3.4–3.5 của bản thảo C đang trình bày kết quả đã biết từ 13 năm trước
như phát hiện mới. **Phải viết lại khung.**

### Phần còn sống sau khi đọc tài liệu

1. **Bất đối xứng khi truyền qua well-stirred.** Cả ba bài dừng ở độ chính xác
   dự đoán. Không bài nào hỏi sai số `fu` đáng giá bao nhiêu trong `CL_h`, và
   không bài nào nhận ra bão hoà chỉ giảm nhẹ cho chất gắn *yếu* — nên vùng sai
   số lớn nhất là vùng không được giảm nhẹ. Chưa thấy ở đâu.
2. **Benchmark tái tạo lại vấn đề đã giải.** TDC `PPBR_AZ` vẫn phát hành thang
   phần trăm và xếp hạng bằng MAE trên phần trăm, tám năm sau Zhu 2013.
3. **Giải thích bằng định giá hậu quả** thay cho lý do thống kê. Zhu và Watanabe
   biện minh bằng "phân bố gần chuẩn hơn"; số học gradient (111×) nói *vì sao*
   điều đó quan trọng. Là giải thích tốt hơn, không phải phát hiện mới.
4. **Sàn MAE của benchmark** (mục 3 dưới) và **gộp loài** (mục 2) — cả hai mới.

## 2. `PPBR_AZ` gộp năm loài, 42,4% không phải người

Tôi tìm ra cái này tình cờ: metadata của TDC ghi `"ppbr_az": 1797` nhưng
`admet_group` trả về 2790 hàng. Đáng lẽ phải kiểm từ đầu.

Deposition `CHEMBL3301361` chứa **năm** assay gắn protein huyết tương:

| ChEMBL | Loài | Bản ghi |
|---|---|---|
| CHEMBL3301365 | người | 1614 |
| CHEMBL3301366 | chuột cống | 717 |
| CHEMBL3301367 | chó | 244 |
| CHEMBL3301368 | chuột nhắt | 162 |
| CHEMBL3301369 | chuột lang | 91 |
| | **tổng** | **2828 trên 1797 phân tử** |

`PPBR_AZ`: 2790 hàng trên **1797 phân tử**. Số phân tử khớp chính xác.

Khớp giá trị theo đúng phương pháp paper A mục 4.1 (97,5% hàng khớp đúng một
assay):

| Loài | Hàng | Tỉ lệ (trong số khớp rõ) |
|---|---|---|
| người | 1567 | 57,6% |
| chuột cống | 680 | 25,0% |
| chó | 233 | 8,6% |
| chuột nhắt | 153 | 5,6% |
| chuột lang | 86 | 3,2% |
| | | **42,4% không phải người** |

Đây là khuyết điểm gộp loài **thứ hai** trong cùng benchmark, trên endpoint thứ
hai, với **năm** loài thay vì hai. Paper A mạnh lên đáng kể: không còn là một sơ
suất ở một endpoint mà là một kiểu chuyển đổi.

### Mệnh đề sai trong paper A, phải sửa trước khi nộp

Paper A mục 4.1 viết: "`Clearance_Microsome_AZ` and `PPBR_AZ` match a single
human assay exactly, at 1 102 and 1 614 records."

- Microsome: benchmark 1102 (881 + 221) = ChEMBL 1102. **Đúng.**
- PPB: benchmark 2790 hàng / 1797 phân tử, ChEMBL người 1614. **Sai.**

Câu đó lẫn "assay ngược dòng có 1614 bản ghi" với "benchmark khớp nó chính xác".

### Và nó giải thích hậu nghiệm cái trần nhiễu đã thu hồi

Paper A mục 4.4 thu hồi một trần nhiễu ước từ độ trải trong-phân-tử, vì các
"nhóm lặp" của hepatocyte là cặp người/chuột cống. PPB có 691 phân tử lặp với
độ trải tới 62,77 điểm — và ngược dòng, assay người có **1614 bản ghi trên 1614
phân tử, không trùng lặp một cái nào**. Nên các nhóm lặp của PPB cũng là khác
loài, không phải đo lặp. Lần thứ hai, cùng một cái bẫy, cùng một dataset. Bài
học của mục 4.4 tổng quát hơn tôi tưởng khi viết nó.

## 3. Sàn MAE không thể vượt của bảng xếp hạng

Test set của `PPBR_AZ` có **559 hàng nhưng 343 phân tử**; 140 phân tử mang nhiều
nhãn mâu thuẫn (trải tới 49,27 điểm). Dự đoán tốt nhất có thể cho một phân tử
dưới MAE là trung vị nhãn của nó, nên tồn tại sàn:

**Sàn MAE = 1,964 điểm phần trăm. Sàn RMSE = 4,478.**

Đối chiếu: đầu bảng xếp hạng ở 7,44 (bản sao bên thứ ba, chụp 2026-05-28);
`logit` của tôi cho 7,433; `percent` cho 7,812. Các mô hình đang tranh nhau
**0,379 điểm** trong khi **1,964 điểm** không mô hình nào bỏ được — và sàn đó
chính là bất đồng giữa các loài, không phải nhiễu assay.

Không trùng lặp nào vượt ranh giới train/test (0 phân tử), nên mệnh đề "zero
overlap" của paper A mục 4.3 vẫn đứng và không có rò rỉ.

## 4. Tách hiệu ứng thang khỏi hiệu ứng loài

`src/experiments/human_only_scale.py`. Huấn luyện trên bộ gộp như phát hành so
với chỉ assay người (1614 bản ghi, một bản ghi một phân tử), đánh giá trên 128
chất ngoài miền và trên phần người của test set. 25 seed, ghép cặp.

Trung vị sai số `fu`, nhóm gắn >99%, trên bộ ngoài miền:

| Huấn luyện | `percent` | `log_fu` | `logit` |
|---|---|---|---|
| Gộp loài (như phát hành) | 8,710 | 2,351 | 2,458 |
| Chỉ người | 8,779 | **2,158** | **2,154** |

**Dự đoán 1 đúng — hiệu ứng thang sống sót:** −6,62 lần ở bộ chỉ người so với
−6,36 ở bộ gộp, 25/25 seed cả hai. Luận điểm thang không phải sản phẩm phụ của
gộp loài.

**Dự đoán 2 của tôi sai.** Tôi dự đoán bỏ dữ liệu khác loài sẽ giảm sai số `fu`
ở mọi thang. Ở `percent` nó **không đổi gì**: +0,070, p=0,85, 13/25. Ở hai thang
biến đổi chỉ giảm 0,19–0,30.

Nên: **hai khuyết điểm chưa công bố trong cùng một endpoint, và cái lớn hơn
không phải cái dễ thấy.** Gộp năm loài là lỗi tài liệu hiển nhiên và gần như
không tốn gì về độ chính xác `fu`; sai thang mục tiêu không ai gọi là khuyết
điểm và tốn một hệ số bốn. Đây là kết quả tốt nhất của phiên, và nó đến từ một
dự đoán sai.

(Lần thứ ba trong dự án tôi dự đoán sai về một đánh đổi. notes/09 Tobit, notes/21
thang, và đây. Ba lần đều theo chiều "tôi tưởng hai thứ liên quan hơn thực tế".)

## 5. Việc tiếp theo

- Sửa mệnh đề 1614 trong paper A và thêm phát hiện gộp năm loài PPB. Paper A
  đang chờ nộp nên việc này ưu tiên cao nhất.
- Viết lại khung paper C quanh mục 1 (bất đối xứng), 2 (gộp loài), 3 (sàn MAE),
  4 (tách hai hiệu ứng), hạ phần thang xuống thành tái lập + phê bình benchmark
  và trích Zhu/Watanabe/Han đúng mức.
- Cập nhật issue #405 gửi TDC: nó hiện chỉ nêu hepatocyte và kiểm duyệt. Phải
  thêm PPB gộp năm loài và chênh lệch 2790/1797/1614.
- Ba hình cho C, và quyết định việc "so với mô hình đã công bố" — giờ câu hỏi
  đúng là bảng xếp hạng có tiếp nhận thang biến đổi hay không, chứ không phải
  chạy mô hình của người khác.
