# Finding 13 — Phương sai giao thức, đo lại bằng đúng giao thức của bảng xếp hạng

Ngày: 2026-10-17. Script: `src/eval/leaderboard_protocol.py`.
Kết quả: `results/leaderboard_protocol_*.csv`. **Thay thế các con số ở `notes/02`.**

## Vì sao phải đo lại

`notes/02` dùng `ADME(...).get_split(method='scaffold', seed=s)`, chia lại **toàn bộ** dữ
liệu kể cả test set trên mỗi seed. Giao thức thật của `admet_group` có **test set cố
định** (243 phân tử cho `cl_hep`, hash không đổi theo seed); seed chỉ chia lại train/valid
(848/122). Đã tự kiểm chứng bằng hash.

Hai giao thức đo hai thứ khác nhau, nên đo cả hai trên cùng mô hình và cùng đặc trưng.

## Kết quả, 25 seed mỗi bên

| Endpoint | Chỉ số | **Giao thức bảng** (test cố định) | | **Chia lại cả test** | | tỷ lệ SD |
|---|---|---|---|---|---|---|
| | | TB | **SD** | TB | SD | |
| `cl_hep` | Spearman | 0,4406 | **0,0201** | 0,4896 | 0,0538 | **2,67×** |
| `cl_mic` | Spearman | 0,5569 | **0,0213** | 0,5773 | 0,0506 | **2,38×** |
| `vdss` | Spearman | 0,6534 | **0,0142** | 0,6857 | 0,0451 | **3,18×** |
| `half_life` | Spearman | 0,4616 | **0,0247** | 0,4520 | 0,0865 | **3,51×** |
| `ppbr` | MAE | 7,9231 | **0,1848** | 8,6174 | 0,7485 | **4,05×** |

SD dưới giao thức thật là **0,014–0,025**, cùng bậc với 0,005–0,016 mà Koleiev et al. báo
(khác biệt hợp lý vì khác mô hình). Con số 0,045–0,087 ở `notes/02` lớn hơn **2,4–4,1
lần**, và chênh đó **được giải thích**: nó gồm thêm phương sai của việc *chọn phân tử nào
để chấm điểm*.

## Phát biểu đã sửa

### Cái gì SAI ở notes/02

*"Cùng một mô hình chạy hai lần theo giao thức của bảng chênh nhau ≥0,067 ở 5% số lần."*
→ Con số đó thuộc giao thức **chia lại cả test**, không phải giao thức của bảng.

### Cái gì ĐÚNG, đo bằng giao thức thật

Cùng một mô hình, hai lần chạy 5 seed rời nhau, **test set cố định**:

| Endpoint | chênh lệch giả, phân vị 95 |
|---|---|
| `cl_hep` | **0,0247** |
| `cl_mic` | **0,0260** |
| `vdss` | 0,0174 |
| `half_life` | **0,0302** |
| `ppbr` (MAE) | 0,2272 |

Nhỏ hơn con số cũ nhiều, **nhưng vẫn đủ lớn để nuốt phần lớn khoảng cách trên bảng**:

- Top-5 của `Clearance_Microsome` chênh nhau **0,005**, trong khi chênh lệch giả của chính
  endpoint này là **0,026** — gấp **năm lần** toàn bộ khoảng cách top-5.
- Bảng 3 của Koleiev et al. cho thấy `VDss` và `Half_Life` có **hai mô hình đổi chỗ theo
  hai chiều ngược nhau** trên chênh lệch 0,001–0,009 — nằm gọn trong 0,017–0,030.
- Thí nghiệm overfit cố ý của họ đưa `Clearance_Microsome` từ hạng 14 lên hạng 1 với mức
  tăng **+0,077**, và các số đó **chạy một lần, không thanh sai số**.

### Và một phát biểu mới, mạnh hơn cái cũ

> **Dấu ± mà bảng công bố chỉ gồm phương sai phân hoạch train/valid và tính ngẫu nhiên của
> bộ học, vì test set cố định. Lấy mẫu lại chính các phân tử được chấm điểm làm phương sai
> tăng 2,4–4,1 lần. Nên dấu ± trên bảng không nói gì về việc con số đó sẽ dịch bao nhiêu
> trên một mẫu phân tử khác — tức đúng câu hỏi mà người đọc bảng quan tâm.**

Đây là phát hiện tốt hơn cái ban đầu, vì nó **giải thích được** tại sao số của mình và số
của Koleiev lệch nhau 4–8 lần, thay vì mâu thuẫn với họ.

## Lưu ý khi so sánh với mốc đầu bảng

Dưới giao thức thật, baseline của mình huấn luyện trên phân hoạch `train` của từng seed
(848 phân tử), nên trung bình 0,4406 thấp hơn con số 0,4896 ở `notes/02` vốn huấn luyện
trên train+valid (970). So với mốc tái lập được MapLight+GNN **0,500 ± 0,012**
(`notes/12`), baseline của mình kém khoảng **0,06**, tức ~3 SD. **Baseline thật sự yếu
hơn** — không được đọc phát hiện phương sai thành ngang bằng.

## Giới hạn

1. Hai giao thức khác nhau cả ở lượng dữ liệu huấn luyện (848 so với 970), nên không tách
   riêng được phần đóng góp của việc chia lại test khỏi phần đóng góp của ít dữ liệu hơn.
   Muốn tách sạch thì phải cố định cỡ tập huấn luyện ở cả hai nhánh.
2. Chỉ một mô hình (LightGBM trên descriptor). Phương sai có thể khác với mô hình khác.
3. Koleiev dùng 5 seed; mình dùng 25. SD của họ do đó kém tin cậy hơn, đúng như `notes/09`
   đã chứng minh cho chính trường hợp 5 seed.
