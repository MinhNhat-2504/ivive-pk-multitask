# Finding 01 — Độ phủ nhãn giữa 5 endpoint PK

Ngày đo: 2026-10-17. Script: `src/data/load_tdc.py`. Dữ liệu: TDC ADMET, PyTDC 0.4.1,
SMILES đã canonical hoá bằng RDKit rồi gộp trùng (trung bình các lần đo lặp).

## Số liệu

Kích thước từng endpoint sau khi canonical hoá + gộp trùng:

| key | TDC dataset | n |
|---|---|---|
| `cl_hep` | Clearance_Hepatocyte_AZ | 1 020 |
| `cl_mic` | Clearance_Microsome_AZ | 1 102 |
| `ppbr` | PPBR_AZ | 1 614 |
| `vdss` | VDss_Lombardo | 1 111 |
| `half_life` | Half_Life_Obach | 665 |

Hợp của cả năm: **3 898 phân tử duy nhất**.

Đồng xuất hiện nhãn theo cặp:

|  | cl_hep | cl_mic | ppbr | vdss | half_life |
|---|---|---|---|---|---|
| **cl_hep** | 1020 | 682 | 521 | **25** | **93** |
| **cl_mic** | 682 | 1102 | 541 | **16** | **54** |
| **ppbr** | 521 | 541 | 1614 | **38** | 119 |
| **vdss** | 25 | 16 | 38 | 1111 | 100 |
| **half_life** | 93 | 54 | 119 | 100 | 665 |

Số phân tử theo số endpoint có nhãn: 1 nhãn → 2 791; 2 → 660; 3 → 395; 4 → 44; **5 → 8**.

## Diễn giải: dữ liệu chia làm hai nhánh gần như rời nhau

- **Nhánh in vitro** (`cl_hep`, `cl_mic`, `ppbr`) nối với nhau rất tốt — cùng nguồn
  AstraZeneca. cl_hep ∩ cl_mic = 682; mỗi cái ∩ ppbr ≈ 520–540.
- **Nhánh in vivo** (`vdss`, `half_life`) nối với nhau tạm được (100) — Lombardo và Obach.
- **Cầu nối giữa hai nhánh chỉ 16–38 phân tử.** Đây đúng là chỗ IVIVE sống.

Hệ quả trực tiếp:

| Ràng buộc | Cần | Số phân tử |
|---|---|---|
| Well-stirred: CL_hep + f_u → CL_int,in vivo | cl_hep + ppbr | **521** |
| Well-stirred: CL_mic + f_u → CL_int,in vivo | cl_mic + ppbr | **541** |
| t½ = ln2 · VDss / CL (khép vòng) | half_life + vdss + cl_hep | **20** |
| Toàn chuỗi, đủ 5 nhãn | cả 5 | **8** |

## Kết luận — cách làm hiển nhiên đã chết, và đó là cơ hội

Cách làm hiển nhiên là đặt ràng buộc IVIVE thành **một số hạng loss có giám sát** trên
các phân tử có đủ nhãn. Với 20 phân tử thì không làm được. Nếu cứ làm, kết quả sẽ nằm
gọn trong sai số và đúng y như lời phê bình đã cảnh báo.

Nhưng ràng buộc vật lý **không cần nhãn mới áp được**. Với bất kỳ phân tử nào, mô hình
đa nhiệm vẫn xuất đủ 5 đầu ra, và quan hệ well-stirred + t½ = ln2·Vd/CL phải đúng
**giữa các giá trị dự đoán với nhau**, bất kể phân tử đó có nhãn hay không. Vậy ràng
buộc áp được trên toàn bộ 3 898 phân tử — và mở rộng được sang phân tử không nhãn lấy
từ ChEMBL.

Điều này đảo ngược luận điểm của bài:

> Nhãn PK thưa và phân mảnh theo nguồn đo là **lý do** phải nhúng dược động học vào mô
> hình, chứ không phải trở ngại. Ràng buộc cơ chế là thứ truyền được thông tin qua chỗ
> nhãn không bắc cầu.

Và **20 phân tử đủ 3 nhãn khép vòng trở thành tập kiểm định ràng buộc**, không phải dữ
liệu huấn luyện — đúng vai trò của nó.

So với công trình đã có: PEMAL dùng hằng số cố định (K₁=16,846; K₂=17,71) trên chuột;
MTGBM không có ràng buộc công thức nào. Chưa ai áp ràng buộc IVIVE kiểu bán giám sát
với phần thuốc tự do tự học.

## Vấn đề kỹ thuật còn mở

Gradient boosting không nhận trực tiếp một loss khớp nối nhiều đầu ra. Hai lối:

1. **Vòng lặp refit có ràng buộc** — fit các head bằng LightGBM, chiếu dự đoán lên mặt
   thoả ràng buộc, refit theo pseudo-target, lặp. Giữ được lợi thế boosting mà PDF chỉ
   ra (boosting thắng GNN trên chính các endpoint này), chạy CPU.
2. **MLP nhỏ trên descriptor** — loss khớp nối khả vi trực tiếp, nhưng bỏ mất lợi thế
   boosting.

Đề xuất: làm (1) là phương pháp chính, (2) là ablation.

## Việc tiếp theo

- [ ] Baseline đơn nhiệm LightGBM trên split scaffold TDC, 5 seed — xác nhận tái lập
      được 0,536 / 0,576 trước khi cải tiến bất cứ thứ gì.
- [ ] Ước lượng trần nhiễu từ các lần đo lặp trong chính dữ liệu thô.
- [ ] Dựng external set theo thời gian từ Drugs@FDA (thuốc duyệt sau 2016).
