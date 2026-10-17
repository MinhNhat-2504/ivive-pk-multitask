# Finding 06 — Ràng buộc IVIVE đúng tới đâu trên dữ liệu thật

Ngày đo: 2026-10-17. Module: `src/models/ivive.py`. Kiểm **trước khi** dựng mô hình ràng
buộc, vì nếu quan hệ không đúng ngay cả với nhãn thật thì mọi thứ xây lên đều vô nghĩa.

Hằng số sinh lý người: hepatocellularity 99×10⁶ tế bào/g gan, protein microsome 32 mg/g
gan, khối lượng gan 25,7 g/kg, lưu lượng máu gan 20,7 mL/phút/kg. Giả định B:P = 1.

## 1. Chuỗi đầy đủ: CLint + PPBR + VDss → t½

| Chuỗi | n | bias log10(quan sát/IVIVE) | sd | sai số bội trung vị | **Spearman [KTC 95%]** |
|---|---|---|---|---|---|
| Hepatocyte, tất cả | 16 | +0,21 | 0,55 | 2,6× | **0,399 [−0,12 – 0,75]** |
| Hepatocyte, bỏ kiểm duyệt | 13 | +0,10 | 0,56 | 2,0× | 0,355 [−0,24 – 0,76] |
| **Microsome, tất cả** | **10** | **−0,38** | 0,39 | 2,3× | **0,891 [0,60 – 0,97]** |
| Microsome, bỏ kiểm duyệt | 6 | −0,42 | 0,50 | 2,8× | 0,829 [0,05 – 0,98] |

### Đọc ra

**Nhánh microsome hoạt động.** Spearman 0,891, KTC không chứa 0. Có thiên lệch hệ thống
−0,38 log10 → IVIVE **ước lượng thiếu độ thanh thải khoảng 2,4 lần**. Đây đúng là hiện
tượng under-prediction kinh điển mà y văn IVIVE đã ghi nhận hàng chục năm cho microsome.
Việc nó hiện ra đúng dấu và đúng độ lớn là **kiểm chứng độc lập rằng module quy đổi của
ta viết đúng**.

**Nhánh hepatocyte không chứng minh được là hoạt động.** Spearman 0,399 với KTC
[−0,12 – 0,75] **chứa 0**. Không kết luận được là nó hỏng — n = 16 thì CI rộng tới mức
gần như không loại trừ được gì — chỉ kết luận được là **chưa có bằng chứng**.

**Kiểm duyệt không phải nguyên nhân.** Bỏ 3 phân tử bị kiểm duyệt đi, Spearman gần như
không đổi (0,399 → 0,355).

## 2. Ràng buộc hepatocyte ↔ microsome (n = 682, mẫu lớn nhất có được)

Cả hai assay đo cùng một đại lượng — độ thanh thải nội tại — bằng hai hệ thống khác nhau.
Sau khi quy đổi sinh lý, chúng phải khớp.

| | n | bias log10(hep/mic) | sd | **Spearman** |
|---|---|---|---|---|
| Tất cả | 682 | +0,58 | 0,50 | **0,642** |
| **Bỏ kiểm duyệt** | **395** | +0,52 | 0,46 | **0,412** |

Hai điều:

1. **Thiên lệch hệ thống 3,8 lần** (10^0,58): CLint quy đổi từ tế bào gan cao hơn hẳn từ
   microsome. Nhất quán, đo trên 682 phân tử — đây là một hằng số hiệu chỉnh **học được**.
2. **Spearman tụt từ 0,642 xuống 0,412 khi bỏ giá trị kiểm duyệt.** Phần lớn tương quan
   biểu kiến đến từ việc hai assay cùng báo "<3" hoặc cùng báo ">150". Mức đồng thuận
   thật giữa hai phép đo cùng một đại lượng chỉ là **0,41**.

Điểm thứ hai đáng chú ý: **mô hình của ta dự đoán `cl_hep` (Spearman 0,49) tốt hơn mức
mà assay microsome dự đoán được `cl_hep` (0,41).** Hai con số không cùng tập kiểm tra nên
không so trực tiếp được, nhưng cùng cỡ — và nó cho thấy phần "tín hiệu" còn lại của
endpoint này nhỏ đến mức nào.

## 3. Hệ quả cho thiết kế — ràng buộc có thật nhưng YẾU

Sai số của chính ràng buộc là **2–3 lần** (trung vị fold-error 2,0–2,8×). Nhiễu phép đo
ước lượng ở `notes/04` là **2,2 lần** (sd_within 0,798 trên thang log tự nhiên ở nhóm
không kiểm duyệt → e^0,798 ≈ 2,2). **Hai con số cùng cỡ.**

Nghĩa là ràng buộc IVIVE không mang thêm thông tin vượt quá mức nhiễu của chính dữ liệu.
Không được kỳ vọng nó cải thiện điểm trên bảng. Cái nó có thể làm:

1. **Truyền thông tin qua chỗ nhãn không bắc cầu** — đúng vấn đề ở `notes/01`.
2. **Ràng buộc dự đoán vào vùng khả dĩ về mặt sinh lý**, nhất là ngoài miền huấn luyện.
3. **Giảm nhu cầu nhãn** — cùng độ chính xác với ít nhãn hơn.

Cả ba đều **đo được**, và không cái nào là "điểm cao hơn trên TDC".

## 4. Thiết kế đã chốt cho mô hình ràng buộc

- Ràng buộc chính đặt trên **cặp hepatocyte ↔ microsome** (n = 682, mẫu lớn nhất), không
  phải chuỗi đầy đủ (n = 16).
- Chuỗi đầy đủ dùng làm **ràng buộc phụ** và **tập kiểm định**, không phải nguồn tín hiệu.
- **Hệ số hiệu chỉnh học được cho từng hệ thống** (chặn + độ dốc trên thang log) thay cho
  hằng số cố định — đây chính là điểm khác PEMAL, và dữ liệu vừa cho thấy hai hằng số cần
  học: 3,8× (hep/mic) và 2,4× (under-prediction của microsome).
- Trọng số của ràng buộc phải đặt theo **độ tin cậy đã đo được**, không chọn tuỳ ý.
- Bỏ giá trị kiểm duyệt khi tính dư của ràng buộc — chúng làm tương quan giả tăng.

## 5. Giới hạn — phải viết vào bài

1. n = 16 và n = 10 cho chuỗi đầy đủ. Mọi phát biểu về nó đều yếu.
2. Giả định B:P = 1 sai với chất phân bố vào hồng cầu.
3. Bỏ qua hiệu chỉnh gắn kết trong ống nghiệm (fu,inc) — sẽ đưa vào như tham số học được.
4. Bỏ qua thải trừ ngoài gan và qua thận; với chất thải chủ yếu qua thận thì chuỗi sai
   về nguyên tắc, không phải do nhiễu.
5. Hằng số sinh lý là giá trị trung bình quần thể.

## Việc tiếp theo

- [ ] Dựng mô hình đa nhiệm có ràng buộc hep↔mic + hệ số hiệu chỉnh học được.
- [ ] So ghép cặp theo seed với baseline đơn nhiệm.
- [ ] Thí nghiệm giảm nhãn: huấn luyện với 25/50/75% nhãn, xem ràng buộc bù được bao nhiêu.
