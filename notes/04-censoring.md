# Finding 04 — Assay AZ bị kiểm duyệt hai đầu

> ## ⛔ PHẦN TRẦN NHIỄU ĐÃ RÚT LẠI — XEM `notes/15`
>
> Mục 2 và 3 của trang này "sửa" trần nhiễu từ 0,541 lên ≥0,764. **Cả hai con số đều
> không dùng được**: các nhóm đo lặp mà cả hai dựa vào là **cặp người–chuột**, vì
> `Clearance_Hepatocyte_AZ` của TDC gộp hai loài. Không còn ước lượng trần nhiễu nào cho
> `cl_hep`.
>
> **Phần kiểm duyệt (mục 1) vẫn đúng** và giờ có nguồn trích dẫn trực tiếp từ mô tả assay
> của ChEMBL.

Ngày đo: 2026-10-17. Script: `src/eval/censoring.py`. Kết quả: `results/censoring*.csv`.

Manh mối ban đầu đến từ báo cáo tra cứu external set: assay thải trừ nội tại của
AstraZeneca có dải khai báo "<3 đến >150 µL/min/10⁶ tế bào". Kiểm tra trên dữ liệu thật
thì đúng, và hệ quả lớn hơn dự kiến.

## 1. Một phần ba nhãn không phải phép đo

| Endpoint | bản ghi | tại giới hạn dưới (=3) | tại giới hạn trên (=150) | **tổng bị kiểm duyệt** |
|---|---|---|---|---|
| `cl_hep` | 1 213 | 195 (16,1%) | 137 (11,3%) | **27,4%** |
| `cl_mic` | 1 102 | 287 (26,0%) | 84 (7,6%) | **33,7%** |

min = 3,0 và max = 150,0 ở cả hai endpoint, không có giá trị nào ngoài dải. Đây là cờ
"dưới/trên ngưỡng phát hiện" được ghi thành số, rồi được cả cộng đồng huấn luyện như
biến hồi quy liên tục.

Phần thông tin thứ hạng bị mất do hai khối giá trị bằng nhau: 2,9% số cặp ở `cl_hep`,
7,3% ở `cl_mic`. Nhỏ hơn mình đoán, nhưng không bỏ qua được.

## 2. Điều này làm hỏng ước lượng trần nhiễu ở `notes/03` — và mình đã báo cáo sai

Hai sự kiện trong dữ liệu:

- `cl_hep` có **0 bản ghi trùng lặp chính xác cặp (SMILES, giá trị)**.
- Trong 193 nhóm đo lặp, có **0 nhóm mà mọi giá trị đều bị kiểm duyệt**.

Hai điều đó cùng chỉ về một kết luận: **TDC (hoặc nguồn thượng nguồn) đã gộp bỏ các cặp
đo lặp trùng khớp.** Nên 193 "nhóm đo lặp" mình dùng ở `notes/03` không phải mẫu ngẫu
nhiên của các phép đo lặp — chúng là **đúng những cặp bất đồng**. Các phân tử đo hai lần
ra kết quả giống nhau đã biến mất khỏi tầm nhìn.

Thêm nữa, 87/193 nhóm chứa ít nhất một giá trị bị kiểm duyệt. Khi một lần đo là 3,0
(cờ "<3") còn lần kia là giá trị thật, độ tản trông rất lớn — nhưng đó là **hiện vật của
kiểm duyệt, không phải nhiễu phép đo**.

### Độ nhạy

| Tập nhóm đo lặp | n nhóm | sd_within | **ICC** | **trần Spearman** |
|---|---|---|---|---|
| Tất cả nhóm bất đồng (= `notes/03`) | 193 | 1,024 | 0,294 | **0,541** |
| **Chỉ nhóm không có giá trị kiểm duyệt** | **106** | **0,798** | **0,571** | **0,764** |

**Loại bỏ hiện vật kiểm duyệt, ICC tăng từ 0,29 lên 0,57 và trần Spearman từ 0,54 lên
0,76.** Và 0,764 **vẫn là cận dưới**, vì vẫn còn điều kiện "chỉ gồm cặp bất đồng".

## 3. Kết luận đã sửa — ngược với điều ghi ở notes/03 lần đầu

| | Spearman |
|---|---|
| Baseline của ta (25 seed) | 0,490 |
| Đầu bảng TDC | 0,536 |
| ~~Trần nhiễu (bản sai)~~ | ~~0,541~~ |
| **Trần nhiễu (cận dưới, đã sửa)** | **≥ 0,764** |

**`cl_hep` chưa bão hoà. Còn dư địa thật, cỡ ≥ 0,22 Spearman từ mức đầu bảng.**

Điều mình nói với người dùng ở lượt trước — "không phân biệt được mô hình đầu bảng với
mức tối đa đạt được" — **là sai**. Nó đúng với con số 0,541, nhưng con số đó là hiện vật
của kiểm duyệt chứ không phải tính chất của phép đo. Mục tiêu 0,536 → 0,60 trong PDF ban
đầu **nằm dưới trần**, hoàn toàn khả thi về mặt dữ liệu.

Bài học phương pháp, đáng viết vào bài: **ước lượng trần nhiễu từ dữ liệu đã bị gộp trùng
và bị kiểm duyệt sẽ cho trần thấp giả tạo** — tức sẽ kết luận nhầm rằng một benchmark đã
chết trong khi nó chưa.

## 4. Những gì ở `notes/02` KHÔNG bị ảnh hưởng

Phát hiện về phương sai giao thức độc lập hoàn toàn với chuyện này và vẫn đứng nguyên:

- sd giữa các seed 0,045–0,087 Spearman.
- Năm seed không ước lượng nổi chính sd đó.
- Cùng một mô hình chạy hai lần theo giao thức chênh nhau ≥ 0,067 (`cl_hep`) ở 5% số lần.

Vậy kết luận tổng hợp giờ là: **endpoint còn dư địa thật, nhưng thước đo hiện tại không
đo nổi dư địa đó.** Đó là một luận điểm mạnh hơn và cũng dễ bảo vệ hơn.

## 5. Mở ra một hướng cải tiến có cơ sở

27–34% nhãn bị kiểm duyệt mà **chưa ai xử lý như dữ liệu kiểm duyệt**. Hồi quy thường
trên nhãn kiểm duyệt là sai mô hình thống kê: nó ép mô hình khớp giá trị 3,0 như thể đó
là phép đo, trong khi sự thật chỉ là "nhỏ hơn 3".

Hướng: **hàm mất mát có nhận thức kiểm duyệt** (kiểu Tobit — với nhãn bị kiểm duyệt dưới,
chỉ phạt khi dự đoán vượt ngưỡng). Đây là cải tiến có cơ sở thống kê, chạy được trên
LightGBM bằng objective tuỳ biến, và **độc lập với ràng buộc IVIVE** nên làm được ablation
sạch.

## 6. Giới hạn còn lại — phải ghi vào bài

1. Trần 0,764 vẫn là cận dưới: vẫn điều kiện trên cặp bất đồng. Không bỏ được điều kiện
   này nếu không có dữ liệu AZ gốc chưa gộp trùng.
2. `var_between` tính trên **toàn bộ** phân tử (gồm cả nhóm bị kiểm duyệt, vốn bị nén về
   hai đầu) trong khi `var_within` chỉ từ nhóm không kiểm duyệt. Hai phần lệch nhau về
   tập hợp; cần kiểm tra lại bằng cách tính `var_between` trên riêng phần không kiểm duyệt.
3. `cl_mic` chưa ước lượng được trần (0 nhóm đo lặp hiện ra).

## Việc tiếp theo

- [ ] Tính lại `var_between` chỉ trên phần không bị kiểm duyệt (giới hạn 2 ở trên).
- [ ] Thử xin dữ liệu AZ gốc chưa gộp trùng từ ChEMBL (CHEMBL3301361) để bỏ điều kiện bất đồng.
- [ ] Baseline có hàm mất mát nhận thức kiểm duyệt, so ghép cặp theo seed với baseline thường.
