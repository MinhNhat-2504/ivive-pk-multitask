# Finding 07 — Hàm mất mát Tobit không giúp gì. Kết quả âm tính có lực thống kê.

Ngày đo: 2026-10-17. Scripts: `src/models/censored.py`, `src/experiments/censored_vs_plain.py`.
Kết quả: `results/censored_vs_plain_*.csv`.

Thiết kế: ghép cặp theo seed, 25 seed, cùng phép chia cho cả hai mô hình. Hai thước đo
(Spearman trên thang đo gốc, MAE trên thang log mà mô hình khớp), hai tập test (toàn bộ,
và chỉ phần không kiểm duyệt).

## Kết quả — tám so sánh, không cái nào có hiệu ứng

| Endpoint | Tập test | Thước đo | plain | tobit | chênh (Tobit tốt hơn = dương) | KTC 95% | thắng | p |
|---|---|---|---|---|---|---|---|---|
| `cl_hep` | toàn bộ | Spearman | 0,4896 | 0,4830 | −0,0067 | [−0,015, +0,001] | 11/25 | 0,16 |
| `cl_hep` | toàn bộ | MAE log | 0,8799 | 0,8825 | −0,0026 | [−0,008, +0,003] | 9/25 | 0,37 |
| `cl_hep` | không k.d. | Spearman | 0,4003 | 0,3968 | −0,0035 | [−0,012, +0,005] | 11/25 | 0,41 |
| `cl_hep` | không k.d. | MAE log | 0,7283 | 0,7319 | −0,0036 | [−0,009, +0,002] | 10/25 | 0,20 |
| `cl_mic` | toàn bộ | Spearman | 0,5773 | 0,5768 | −0,0005 | [−0,007, +0,006] | 9/25 | 0,54 |
| `cl_mic` | toàn bộ | MAE log | 0,8155 | 0,8152 | +0,0003 | [−0,006, +0,006] | 14/25 | 0,96 |
| `cl_mic` | không k.d. | Spearman | 0,4697 | 0,4663 | −0,0034 | [−0,014, +0,008] | 9/25 | 0,28 |
| `cl_mic` | không k.d. | MAE log | 0,7018 | 0,7044 | −0,0027 | [−0,010, +0,004] | 12/25 | 0,54 |

Mọi KTC đều chứa 0. Mọi p đều > 0,16. Số lần thắng 9–14 trên 25, tức quanh mức ngẫu nhiên.

## Giả thuyết của tôi đã sai

Khi thấy Tobit không giúp Spearman, tôi dự đoán: Spearman là thước đo **thứ hạng**, mà
nhãn kiểm duyệt vốn nằm ở hai cực nên OLS đã xếp đúng thứ hạng sẵn; Tobit sửa **giá trị**
nên sẽ thắng trên MAE. Chạy thêm MAE để kiểm.

**Dự đoán sai.** MAE cũng không có hiệu ứng (chênh −0,0036 đến +0,0003). Tobit đơn giản
là không giúp gì cả, trên cả hai loại thước đo.

## Vì sao — giải thích khả dĩ, chưa kiểm chứng

1. **Mức sai lệch do kiểm duyệt nhỏ hơn sai số của chính mô hình.** Nếu giá trị thật của
   một chất "<3" là 1,5 thì log1p(3) = 1,386 so với log1p(1,5) = 0,916 — lệch 0,47. Sai
   số dư của mô hình là ~0,73–0,88 trên cùng thang. Sửa một sai lệch nhỏ hơn nhiễu nền
   thì không đo được.
2. **Bài toán bị giới hạn bởi nhiễu, không phải bởi hàm hợp lý.** ICC ≈ 0,57, mô hình đạt
   Spearman 0,49 so với trần ≥ 0,76. Thứ đang chặn là chất lượng đặc trưng và nhiễu phép
   đo, không phải dạng của hàm mất mát. Đây là giải thích khớp nhất với các kết quả khác
   trong dự án.
3. **Cây quyết định vốn khá bền với chuyện này.** LightGBM khớp hằng số từng khúc; khối
   giá trị tại 3,0 đủ đông để thành cụm lá riêng bất kể hàm mục tiêu nào.

## Phạm vi của kết luận này

**Có quyền nói:** với cách cài đặt Tobit này — σ ước lượng hai bước từ dư không kiểm
duyệt, LightGBM, descriptor RDKit + Morgan — không có cải thiện đo được trên cả hai
endpoint, hai thước đo, dưới thiết kế ghép cặp phân giải được hiệu ứng cỡ 0,01.

**Không có quyền nói:** rằng không phương pháp nhận thức kiểm duyệt nào giúp được. Chưa
thử: ước lượng σ đồng thời với cây, σ phụ thuộc phân tử, hoặc mô hình hoá kiểm duyệt ở
tầng dữ liệu (ví dụ multiple imputation cho giá trị bị kiểm duyệt).

## Giá trị cho bài

Đây vẫn là đóng góp, chỉ khác dấu so với dự kiến:

> 27–34% nhãn thải trừ là cờ ngưỡng phát hiện và không ai mô hình hoá chúng như dữ liệu
> kiểm duyệt. Sửa đúng về mặt thống kê bằng hàm hợp lý Tobit **không mang lại cải thiện
> đo được nào** — củng cố chẩn đoán rằng thứ chặn các endpoint này là nhiễu phép đo và
> chất lượng đặc trưng, không phải dạng mô hình.

Nó cũng tiết kiệm công cho người khác, và nó chứng minh thiết kế ghép cặp của ta **phát
hiện được cái không có** — một kiểm chứng âm tính cần thiết trước khi tin vào kết quả
dương tính của thí nghiệm giảm nhãn.

## Việc tiếp theo

- [ ] Gộp vào bài như một mục ablation, không phải phụ lục.
- [ ] Nếu còn thời gian: thử ước lượng σ đồng thời, để phát biểu phạm vi chặt hơn.
