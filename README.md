# IVIVE-constrained multitask PK prediction

Dự đoán đồng thời 5 endpoint dược động học trên người của TDC ADMET, ràng buộc bằng
phép ngoại suy in vitro – in vivo (IVIVE) với phần thuốc tự do **tự học** thay vì hằng
số cố định.

| key | TDC dataset | vai trò trong chuỗi IVIVE |
|---|---|---|
| `cl_hep` | Clearance_Hepatocyte_AZ | CL nội tại in vitro (tế bào gan) |
| `cl_mic` | Clearance_Microsome_AZ | CL nội tại in vitro (microsome) |
| `ppbr` | PPBR_AZ | → f_u, **biến trung gian** của ràng buộc |
| `vdss` | VDss_Lombardo | thể tích phân bố |
| `half_life` | Half_Life_Obach | t½ — khép vòng |

Ràng buộc: well-stirred `CL_h = Q_h·f_u·CL_int / (Q_h + f_u·CL_int)`, rồi
`t½ = ln2 · VDss / CL`. PPBR không còn là nhiệm vụ rời mà là biến trung gian.

## Trạng thái

Xem `notes/` theo thứ tự, và bản thảo paper ở [paper/draft.md](paper/draft.md). Hai phát hiện đã có:
- [01 — độ phủ nhãn](notes/01-task-overlap.md): chỉ 20 phân tử có đủ nhãn để khép vòng ràng buộc, nên ràng buộc phải áp **bán giám sát** trên dự đoán.
- [02 — phương sai giao thức](notes/02-baseline-variance.md): giao thức 5 seed của bảng TDC tự chế ra "cải tiến" tới 0,067 Spearman mà không đổi gì về mô hình — lớn hơn mục tiêu cải tiến ban đầu của dự án.
- [03 — trần nhiễu](notes/03-noise-ceiling.md): chỉ 29% phương sai nhãn `cl_hep` là tín hiệu thật; không phân biệt được mô hình đầu bảng với mức tối đa đạt được.

## Cài đặt

Cần Python 3.11 (3.12+ chưa có wheel cho stack này).

```bash
py -3.11 -m venv .venv            # Windows
python3.11 -m venv .venv          # macOS
.venv/Scripts/python -m pip install -r requirements.txt   # Windows
.venv/bin/python -m pip install -r requirements.txt       # macOS
```

## Chạy

```bash
.venv/Scripts/python src/data/load_tdc.py   # tải dữ liệu + báo cáo độ phủ nhãn
```

```bash
.venv/Scripts/python -m src.models.baseline            # baseline don nhiem, 5 endpoint x 5 seed
```

```bash
.venv/Scripts/python -m src.eval.noise_ceiling      # tran nhieu tu cac lan do lap
.venv/Scripts/python -m src.eval.protocol_variance  # phuong sai giao thuc
```
