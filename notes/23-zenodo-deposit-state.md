# 23 — Trạng thái deposit Zenodo

Ghi lại vì hai bài phụ thuộc vào nó, và vì tôi đã suýt để trích dẫn trỏ nhầm.

## Hai DOI, không được lẫn

| | DOI | Trỏ tới |
|---|---|---|
| Concept | **`10.5281/zenodo.23462435`** | luôn là phiên bản mới nhất |
| v1 | `10.5281/zenodo.23462436` | cố định bản 19-10-2026, 9 tệp |
| v2 | `10.5281/zenodo.23471468` | cố định bản hiện tại, 12 tệp |

**Cả hai bài trích concept DOI**, kèm chữ "version v2" để nói rõ đã làm việc trên
bản nào. Ban đầu tôi để version DOI của v1 — nếu không sửa, người đọc theo trích
dẫn sẽ tới bản **không có** ba tệp loài mà paper C dựa vào.

Bài học: Zenodo hiển thị version DOI ở chỗ dễ thấy nhất (ô Details), còn concept
DOI nằm trong một dòng chữ nhỏ ở khối Versions. Mặc định sao chép là sai.

## v2 chứa gì

12 tệp, đã đối chiếu checksum với `zenodo_bundle/` sau khi publish — khớp từng
byte. Ba tệp thêm so với v1:

- `ppb_all_species_upstream.csv` — năm assay PPB theo loài, 2828 bản ghi
- `ppb_species_composition.csv` — quy thuộc từng hàng `PPBR_AZ` về assay
- `ppb_metric_floor.csv` — sàn 1,964 điểm dưới chỉ số bảng xếp hạng

`README.md` cũng đổi (mô tả cột cho ba tệp mới, và nói deposit phục vụ hai bài).
Tám tệp còn lại giữ nguyên.

Related works có hai mục: `isDerivedFrom` ChEMBL CHEMBL3301361 và
`isSupplementedBy` repo GitHub.

## Tái dựng

`src/data/make_zenodo_bundle.py` dựng lại toàn bộ `zenodo_bundle/` từ
`data/external/` và `results/`. Nếu phải lên v3, chạy nó rồi tải cả 12 tệp vào
bản nháp rỗng — **không** dùng nút "Import files", vì đường đó đòi nhớ xoá
`README.md` cũ và tôi đã hướng dẫn sai chỗ này một lần.

## Còn phải làm khi nộp bài

Nếu có v3 trước lúc bài được chấp nhận, đổi chữ "version v2" trong sáu tệp bản
thảo bằng `src/figures/note_deposit_version.py`. Concept DOI thì không phải đụng.
