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
- [02 — phương sai baseline](notes/02-baseline-variance.md): độ lệch giữa seed (0,055–0,122 Spearman) lớn hơn mức cải tiến nhắm tới, nên "vượt SOTA trên TDC" không làm luận điểm chính được.

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
