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

Xem `notes/` theo thứ tự, và bản thảo paper ở [paper/draft.md](paper/draft.md).

- [01 — độ phủ nhãn](notes/01-task-overlap.md): chỉ 20 phân tử có đủ nhãn để khép vòng
  ràng buộc, nên ràng buộc phải áp **bán giám sát** trên dự đoán.
- [02 — phương sai giao thức](notes/02-baseline-variance.md): giao thức 5 seed của bảng
  TDC tự chế ra "cải tiến" tới 0,067 Spearman mà không đổi gì về mô hình — lớn hơn mục
  tiêu cải tiến ban đầu của dự án.
- [03 — trần nhiễu](notes/03-noise-ceiling.md): phương pháp ước lượng trần nhiễu.
  **⚠ Con số trong trang này đã bị `notes/04` bác bỏ một phần.**
- [04 — kiểm duyệt assay](notes/04-censoring.md): 27–34% nhãn clearance là cờ ngưỡng phát
  hiện, không phải phép đo. Sửa lại trần nhiễu từ 0,541 lên **≥ 0,764** — tức endpoint
  **chưa bão hoà**, vẫn còn dư địa thật.
- [05 — external set](notes/05-external-set.md): dựng được cho `ppbr` (n≈150–200) và
  `half_life` (n≈140–165); nhánh in vitro **không lấy được từ hồ sơ FDA**.
- [06 — kiểm chứng quan hệ IVIVE](notes/06-ivive-relation-validated.md): nhánh microsome
  hoạt động (Spearman 0,891) với độ lệch đúng như y văn; nhánh hepatocyte chưa chứng minh
  được. Ràng buộc chuyển sang cặp hep↔mic (n=395).
- [07 — Tobit âm tính](notes/07-tobit-negative.md): hàm mất mát nhận thức kiểm duyệt
  **không cải thiện gì**, trên cả Spearman lẫn MAE. Kết quả âm tính có lực thống kê.
- [08 — giảm nhu cầu nhãn](notes/08-label-efficiency.md): **kết quả chính.** Ràng buộc
  cho +0,053 ở 25% nhãn (p=0,006), hệ số thay thế ~1,7×. Ablation cho thấy lợi ích đến từ
  việc **học** phép hiệu chỉnh, không phải từ dữ liệu thêm.
- [09 — độ nhạy](notes/09-sensitivity.md): kết quả trên không phụ thuộc hai tham số chưa
  tinh chỉnh; cả 8 cấu hình đều dương.

## Cài đặt

Cần Python 3.11 (3.12+ chưa có wheel cho stack này).

```bash
py -3.11 -m venv .venv                                     # Windows
python3.11 -m venv .venv                                   # macOS
.venv/Scripts/python -m pip install -r requirements.txt    # Windows
.venv/bin/python -m pip install -r requirements.txt        # macOS
```

## Chạy

Thay `.venv/Scripts/python` bằng `.venv/bin/python` trên macOS.

```bash
.venv/Scripts/python -m src.data.load_tdc           # tải dữ liệu + báo cáo độ phủ nhãn
.venv/Scripts/python -m src.models.baseline         # baseline đơn nhiệm (mặc định 5 seed)
.venv/Scripts/python -m src.models.baseline --seeds $(seq 1 25) --tag baseline_lgbm_25
.venv/Scripts/python -m src.eval.noise_ceiling      # trần nhiễu từ các lần đo lặp
.venv/Scripts/python -m src.eval.protocol_variance  # phương sai giao thức (cần bản 25 seed)
.venv/Scripts/python -m src.eval.censoring          # kiểm duyệt assay + độ nhạy của trần nhiễu
.venv/Scripts/python -m src.experiments.censored_vs_plain   # Tobit vs OLS, ghép cặp 25 seed
.venv/Scripts/python -m src.experiments.label_efficiency    # kết quả chính + ablation
.venv/Scripts/python -m src.experiments.sensitivity         # quét pseudo_weight / n_rounds
```
