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

Xem `notes/` theo thứ tự. Bản thảo không nằm trong repo; xem mục Giấy phép
ở cuối.

- [01 — độ phủ nhãn](notes/01-task-overlap.md): chỉ 20 phân tử có đủ nhãn để khép vòng
  ràng buộc, nên ràng buộc phải áp **bán giám sát** trên dự đoán.
- [02 — phương sai giao thức](notes/02-baseline-variance.md): giao thức 5 seed của bảng
  TDC tự chế ra "cải tiến" tới 0,067 Spearman mà không đổi gì về mô hình — lớn hơn mục
  tiêu cải tiến ban đầu của dự án.
- [03 — trần nhiễu](notes/03-noise-ceiling.md): phương pháp ước lượng trần nhiễu.
  **⚠ Toàn bộ hướng này đã bị thu hồi.** Xem `notes/12`: các "nhóm lặp" dùng để
  ước lượng hoá ra là cặp người/chuột, nên cả con số gốc 0,541 lẫn con số "đã
  sửa" 0,764 đều không có cơ sở. Giữ lại hai trang này vì bản thảo mô tả chính
  lần sai đó.
- [04 — kiểm duyệt assay](notes/04-censoring.md): 27–34% nhãn clearance là cờ ngưỡng phát
  hiện, không phải phép đo. Phần kiểm duyệt đứng vững; phần hiệu chỉnh trần nhiễu
  trong trang này thì không, xem trên.
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
- [10 — ngoài miền](notes/10-ood-generalisation.md): split scaffold của TDC **không phải**
  phép kiểm ngoại suy mạnh. Trên split cụm với kiểm soát rò rỉ chặt, luận điểm "học
  coupling thắng cố định coupling" **vững và mạnh lên** (p=0,0001); luận điểm "ràng buộc
  thắng đơn nhiệm" thì **yếu đi nhiều**.

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
.venv/Scripts/python -m src.experiments.ood_generalisation --strict --tag ood_strict

# Tập external theo thời gian (FDA)
.venv/Scripts/python -m src.data.fda_nme            # buoc 1: NME tu 2017
.venv/Scripts/python -m src.data.fda_labels         # buoc 2-3: nhan + cau truc
.venv/Scripts/python -m src.data.external_overlap   # loc + kiem trung lap voi TDC
```

## Giấy phép

Mã nguồn: MIT, xem [LICENSE](LICENSE).

Dữ liệu trong `data/` giữ giấy phép riêng của nguồn. Văn bản nhãn thuốc FDA là
tác phẩm của chính phủ Hoa Kỳ; bản gửi AstraZeneca trên ChEMBL (document
`CHEMBL3301361`) là CC BY-SA 3.0. Bộ dữ liệu đã kiểm chứng và phần chú giải
xuất xứ được lưu trữ riêng tại Zenodo dưới CC BY 4.0, concept DOI
[10.5281/zenodo.23462435](https://doi.org/10.5281/zenodo.23462435).

## Bản thảo

Các bản thảo không nằm trong repo này (thư mục `paper/` bị `.gitignore` loại
trừ). Dữ liệu và mã để tái lập mọi con số trong đó thì có.
