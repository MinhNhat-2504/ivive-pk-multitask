# Finding 03 — Trần nhiễu của phép đo: đầu bảng TDC đã chạm trần

Ngày đo: 2026-10-17. Script: `src/eval/noise_ceiling.py`, kết quả ở
`results/noise_ceiling.csv`.

## Cách đo

TDC lưu nhiều bản ghi cho cùng một phân tử trong cùng một endpoint. Sau khi canonical
hoá SMILES, các bản ghi đó gộp thành nhóm đo lặp, và **độ tản trong nhóm chính là độ lặp
lại của phép đo**.

Phân rã phương sai một chiều trên thang mô hình hoá (log1p với clearance/VDss/t½):

- `sd_within` — độ lệch chuẩn gộp của các lần đo quanh trung bình của chính phân tử đó
- `sd_between` — độ lệch chuẩn giữa các phân tử, **đã trừ** thành phần within
- `ICC = var_between / (var_between + var_within)` — tỷ lệ phương sai nhãn là tín hiệu thật

Trần Pearson = √ICC (hiệu chỉnh suy giảm cổ điển). Trần Spearman không có công thức
đóng nên mô phỏng: lấy mô hình đoán đúng tuyệt đối, tương quan với một lần đo nhiễu.

> **Lỗi đã sửa trong lúc làm:** bản đầu dùng thẳng trung bình quan sát làm "truth". Phần
> lớn phân tử chỉ đo một lần nên "trung bình" của chúng chính là một lần đo nhiễu, và độ
> tản của các trung bình đó bị thổi lên bởi đúng thứ nhiễu đang cần đo. Kết quả là trần
> Spearman ra 0,774 trong khi trần Pearson ra 0,542 — vô lý. Sửa bằng cách co các trung
> bình về trung bình chung cho tới khi phương sai khớp `sd_between`. Hai ước lượng giờ
> khớp nhau (0,541 / 0,542).

## Cấu trúc đo lặp trong dữ liệu thô

| Endpoint | bản ghi | phân tử | nhóm lặp | trùng y hệt giá trị | **lặp thật** |
|---|---|---|---|---|---|
| `cl_hep` | 1 213 | 1 020 | 193 | 0 | **193** |
| `cl_mic` | 1 102 | 1 102 | 0 | 0 | **0** |
| `ppbr` | 1 614 | 1 614 | 0 | 0 | **0** |
| `vdss` | 1 130 | 1 111 | 19 | 2 | **17** |
| `half_life` | 667 | 665 | 2 | 0 | **2** |

Chỉ `cl_hep` có đủ số nhóm để ước lượng tin cậy. Không nhóm nào trùng y hệt giá trị ở
`cl_hep` — tức đây là đo lặp thật, không phải bản ghi nhân đôi.

## Kết quả

| Endpoint | nhóm dùng | sd_within | sd_between | **ICC** [KTC 95%] | **trần Spearman** [KTC 95%] |
|---|---|---|---|---|---|
| `cl_hep` | 193 | 1,024 | 0,660 | **0,294** [0,165 – 0,420] | **0,541** [0,400 – 0,652] |
| `vdss` | 17 | 0,196 | 0,859 | 0,951 [0,911 – 0,979] | 0,936 [0,898 – 0,968] |

KTC từ bootstrap 500 lần, **lấy mẫu lại theo nhóm đo lặp** chứ không theo từng phép đo —
nhóm mới là đơn vị mang thông tin về độ lặp lại, và lấy mẫu trong nhóm cố định sẽ giấu
mất bất định do chỉ có ít nhóm.

## Điều này có nghĩa gì

### 1. Chỉ 29% phương sai nhãn thải trừ tế bào gan là tín hiệu thật

ICC = 0,294. **71% biến thiên trong nhãn `cl_hep` là nhiễu phép đo**, không phải khác
biệt giữa các phân tử.

`sd_within` = 1,024 trên thang log1p. Với giá trị lớn, log1p ≈ ln, nên một độ lệch chuẩn
tương ứng hệ số e^1,024 ≈ **2,8 lần**. Con số này khớp với hiểu biết dược học sẵn có:
thải trừ nội tại đo trên tế bào gan vốn lệch 2–3 lần giữa các lần thí nghiệm. Đây là
kiểm chứng độc lập cho thấy phép phân rã không bị hỏng.

### 2. Không phân biệt được mô hình đầu bảng với mức tối đa đạt được

| | Spearman |
|---|---|
| Baseline LightGBM của ta | 0,502 |
| **Đầu bảng TDC** | **0,536** |
| **Trần nhiễu ước lượng** | **0,541** [KTC 95%: 0,400 – 0,652] |

**Phải phát biểu cẩn thận.** Ước lượng điểm của trần là 0,541, chỉ trên đầu bảng 0,005 —
nhưng khoảng tin cậy rộng [0,400 – 0,652], và cả 0,536 lẫn mục tiêu 0,60 đều nằm **trong**
khoảng đó. Vậy câu đúng là:

> Không phân biệt được mô hình đứng đầu bảng với mức tối đa mà dữ liệu cho phép đạt tới.

**Không** được viết "đầu bảng đã chạm trần" hay "0,60 nằm trên trần" — dữ liệu chưa đủ
để khẳng định. Nhưng câu phát biểu đúng vẫn đủ mạnh, vì nó đảo trách nhiệm chứng minh:
ai công bố cải tiến trên endpoint này phải chứng minh mức cải tiến đó nằm trên trần
nhiễu, chứ không chỉ trên mô hình trước đó.

Thêm nữa, ước lượng điểm của trần (0,541) **thấp hơn** mục tiêu 0,60 đặt ra trong PDF.
Để đạt 0,60 một cách thật sự thì trần phải nằm ở nửa trên của khoảng tin cậy. Đó là một
canh bạc, không phải một kế hoạch.

Độ rộng của khoảng này cũng là lý do phải làm tiếp việc ở mục "Việc tiếp theo": 193 nhóm
chưa đủ để thắt khoảng lại.

### 3. Không có rò rỉ phân tử trong split — kết quả phủ định, nhưng cần ghi

Kiểm tra giao giữa train+valid và test theo SMILES canonical, seed 1–3, trên `cl_hep`,
`vdss`, `half_life`: **giao = 0 ở mọi trường hợp**. Split scaffold của TDC giữ các bản
ghi trùng cùng một phía. Vậy phương sai ở `notes/02` và trần nhiễu ở đây **không phải do
rò rỉ** — chúng là tính chất của dữ liệu.

## Giới hạn của kết luận này — phải viết vào bài

1. 193 nhóm lặp là một tập con; độ nhiễu của chúng chưa chắc đại diện cho cả endpoint.
2. ICC giả định toàn bộ tản trong nhóm là nhiễu. Nếu TDC gộp dữ liệu từ các điều kiện
   thí nghiệm khác nhau thì một phần là hiệu ứng điều kiện. Nhưng với mô hình bị chấm
   điểm trên một nhãn đơn lẻ thì phần đó vẫn không khử được — trần vẫn đúng.
3. `vdss` ICC = 0,951 dựa trên **17 nhóm** — quá ít, không được trích dẫn như một con số
   chắc chắn. Chỉ nên nói "không có bằng chứng VDss bị nhiễu nặng như clearance".
4. `cl_mic`, `ppbr`, `half_life` **không có lần đo lặp nào** → không ước lượng được từ
   trong dữ liệu. Phải lấy ước lượng độ lặp lại từ y văn cho ba endpoint này.

## Hệ quả cho bài

Đây giờ là **kết quả mạnh nhất đang có**, mạnh hơn cả mô hình có ràng buộc:

> Trên Clearance_Hepatocyte_AZ, chỉ 29% phương sai nhãn là tín hiệu thật, và không phân
> biệt được mô hình đứng đầu bảng với mức tối đa dữ liệu cho phép. Cộng đồng đang tối ưu
> một đại lượng mà phần lớn biến thiên của nó là nhiễu phép đo.

Và nó định lại mục tiêu của ràng buộc IVIVE: không phải để **vượt** trần trên endpoint
riêng lẻ — chuyện đó gần như bất khả — mà để **đạt cùng mức chính xác với ít nhãn hơn,
nhất quán cơ chế hơn, và khái quát ra miền ngoài tốt hơn**. Đó mới là thứ có chỗ cải thiện.

## Việc tiếp theo

- [ ] Lấy ước lượng độ lặp lại từ y văn cho `cl_mic`, `ppbr`, `half_life`.
- [ ] Khoảng tin cậy bootstrap cho ICC của `cl_hep`.
- [ ] Kiểm tra độ nhạy: ICC thay đổi thế nào nếu loại các nhóm lệch cực lớn.
- [ ] Đối chiếu baseline 25 seed (`notes/02`) với trần này trên cùng một hình.
