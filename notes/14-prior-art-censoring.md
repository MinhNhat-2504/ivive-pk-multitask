# Finding 14 — Có công trình đi trước về kiểm duyệt. Phải sửa tuyên bố mới lạ.

Ngày: 2026-10-17. Nguồn: đọc toàn văn hai bài audit.

## 1. ⚠ Công trình đi trước đã ghi nhận cờ bất đẳng thức

**Parrondo-Pizarro R, Menestrina L, Garcia-Serna R, Fernández-Torras A, Mestres J.
"Enhancing molecular property prediction through data integration and consistency
assessment." *J Cheminform.* 2025;17:163. doi `10.1186/s13321-025-01103-3`. PMC12574261.**

Trong phần Methods, họ viết:

> *"we excluded referential data included in the TDC clearance dataset. Specifically, we
> retained only intrinsic clearance values from human microsomes **that were not annotated
> with '>' or '<' symbols, as we found these categorical annotations to bias the
> regressor**. This refinement reduced the number of molecules **from 1,102 to 744**."*

1 102 − 744 = **358 phân tử = 32,5% của `Clearance_Microsome_AZ`** mang cờ bất đẳng thức.
Con số này **nằm đúng trong dải 27–34% mình đo**, và 1 102 đúng bằng kích thước
`Clearance_Microsome_AZ`. Cùng một dataset, cùng một hiện tượng.

**Không được tuyên bố "chưa ai nhận ra nhãn clearance AZ bị gắn cờ".** Sai.

> Lý do tra cứu lần đầu không tìm ra: **cộng đồng này không dùng từ "censored".** Grep
> toàn văn cả hai bài cho *censor, detection limit, LLOQ, truncat, saturat, qualifier* →
> **0 kết quả**. Phát hiện chỉ lộ ra khi tìm theo *hiện tượng* (`>`/`<`, "annotated",
> "categorical annotation") và theo **chỗ cỡ mẫu tụt trong Methods**. Ghi lại như một bài
> học tra cứu.

## 2. Nhưng phần của mình vẫn còn, và nhiều hơn — bốn điểm cụ thể

| | Họ | Ta |
|---|---|---|
| Endpoint | **chỉ `Clearance_Microsome_AZ`** (họ chọn microsome vì nhiều chất sống sót hơn sau tiền xử lý) | **cả hai**; `Clearance_Hepatocyte_AZ` họ **không chạm tới** |
| Ngưỡng | không nêu | **3 và 150 µL/phút**, và chỉ ra đống dồn phân bố |
| Tỷ lệ | không phát biểu; chỉ suy ra được bằng cách trừ hai cỡ mẫu trong một câu Methods | **27,4% và 33,7%**, phát biểu thẳng |
| Cách xử lý | **xoá bỏ** — mất một phần ba dữ liệu và làm lệch phân bố còn lại | **hàm hợp lý kiểm duyệt** — dùng được thông tin đó |
| Hệ quả cho trần nhiễu | không bàn | **trần nhiễu ước lượng sai 0,541 thay vì ≥0,764 nếu bỏ qua kiểm duyệt** |

## 3. Và một điểm nữa làm đóng góp mạnh lên — đã tự kiểm chứng

**TDC không ship cờ `>`/`<` gì cả.** Cột `Y` của `Clearance_Microsome_AZ` là `float64`,
không giá trị nào là chuỗi chứa dấu bất đẳng thức. Nên Parrondo-Pizarro phải quay về
ChEMBL mới thấy được cờ.

> Từ TDC một mình, **cách duy nhất phát hiện kiểm duyệt là nhận ra đống dồn ở đúng 3 và
> 150**. Mọi mô hình huấn luyện thẳng trên TDC đều không có cách nào biết.

Đây là một câu phải có trong bài.

## 4. Phát biểu đã sửa — dùng câu này, không dùng câu cũ

> Các chú thích bất đẳng thức trong bộ clearance microsome đã được ghi nhận **một lần,
> thoáng qua**, và xử lý bằng cách xoá bỏ [Parrondo-Pizarro 2025]. Cấu trúc ngưỡng phát
> hiện, độ lớn của nó trên **cả hai** endpoint clearance, việc TDC không phát hành cờ nên
> chỉ phát hiện được qua đống dồn phân bố, và hệ quả của nó lên ước lượng trần nhiễu,
> thì chưa được báo cáo.

Trung thực hơn và mạnh hơn một tuyên bố mới lạ trần trụi — và nó cho một chỗ trích dẫn
để neo vào.

## 5. Schuh et al. — mâu thuẫn chưa giải quyết

**Schuh MG, Daniluk A, Sieber SA. "Auditing widely used biomolecular benchmarks reveals
systematic data inconsistencies." *Chem Sci.* 2026;17(36):17455–17467.
doi `10.1039/d6sc01799a`. PMC13425739.** Công cụ: `benchaudit`.

Họ báo **`Clearance_Hepatocyte_AZ` có tỷ lệ nhiễm cao nhất toàn bộ audit: 4,95%**.

Nhưng mình đo trên chính split `admet_group` chính thức, **ba mức chặt dần**:

| Mức khớp | overlap train_val ↔ test |
|---|---|
| SMILES canonical (giữ lập thể) | **0** |
| Bỏ lập thể | **0** |
| Bỏ muối + khung InChIKey | **0** |

Zero ở cả năm endpoint. **Chưa giải thích được.** Giả thuyết khả dĩ nhất: họ audit một
split khác — `get_split()` ngẫu nhiên sẽ rải 193 **bản ghi** trùng của
`Clearance_Hepatocyte_AZ` (1 213 bản ghi / 1 020 phân tử) sang cả hai phía, điều mà phép
đo của mình không thấy vì mình gộp trùng trước.

**Đã hỏi lại, agent đang chạy tiếp.** Không viết gì về con số 4,95% cho tới khi giải quyết
xong — nếu nó cùng split thì **mình có bug**. Điều đã biết: tiêu chí của họ là *"cross-split
contamination (train/valid→test overlap of cleaned SMILES)"*, tức trùng chính xác chứ
không phải lân cận — nên mâu thuẫn là thật, không phải khác định nghĩa.

Điểm quan trọng không phụ thuộc vào việc đó: tiêu chí "label conflict" của họ là **bất
đồng giữa các SMILES trùng nhau**. Một dataset mà 30% nhãn đều ghim đúng tại 150 µL/phút
sẽ hiện ra là **cực kỳ tự nhất quán** theo tiêu chí đó — không có xung đột nào, vì các
bản trùng khớp nhau hoàn hảo. **Audit của họ mù cấu trúc với kiểm duyệt**, và điều đó nói
được trong bài kèm lý do cụ thể.

## 5b. Bổ sung từ SI — ba điểm

**(i) SI Table S2 của họ có cả hepatocyte: 1 245 → 289 (bỏ 77%).**
**KHÔNG được trích 77% như tỷ lệ kiểm duyệt.** Pipeline của họ còn chạy chuẩn hoá, khử
trùng và gộp InChIKey, và con số khởi điểm 1 245 là ChEMBL thô chứ không phải 1 213 đã
gộp của TDC. Mức giảm đó **bị nhiễu bởi nhiều bước**, không quy được cho riêng bộ lọc cờ.
Nó *gợi ý mạnh* rằng assay hepatocyte bị gắn cờ nhiều — khớp với 27,4% mình đo — nhưng
như công bố thì không dùng làm bằng chứng được.

**(ii) Lệch của TDC clearance là lệch thật, không phải khác miền.** Chênh lệch tuyệt đối
chuẩn hoá trên các chất **dùng chung**: TDC ≈ **0,102–0,126**, trong khi mọi cặp nguồn
khác ≈ **0,001–0,003** — chênh ~50 lần. Tức TDC báo clearance khác hẳn cho chính những
phân tử nó chia sẻ với chuẩn vàng.

**(iii) Trần của TDC clearance thấp, đo độc lập.** Huấn luyện trên TDC → kiểm nơi khác cho
R² **−0,68 đến −1,52** (tệ hơn cả đoán trung bình); TDC đứng riêng được **0,19**. Khớp
hướng với công việc trần nhiễu của mình.

## 5c. ✅ Họ **tán thành thẳng** thiết kế của mình — trích dẫn được

Hai câu trong phần khuyến nghị của họ:

> *"in the case of intrinsic clearance, annotations for the fraction unbound in plasma
> (f_ub) can be used to estimate hepatic metabolic clearance via **in vitro–in vivo
> extrapolation (IVIVE)** methods"*

> *"adopt a **multi-task setting**, where divergent datasets are treated as independent
> tasks rather than a single aggregated dataset… improving both predictive performance
> and generalization"*

Đây là **tán thành đã qua bình duyệt, trên J Cheminform, cho đúng kiến trúc của bài**:
IVIVE dùng phần thuốc tự do, cộng đa nhiệm, thay vì gộp dữ liệu một cách ngây thơ. Đưa
vào Introduction.

Họ còn cho một heuristic chẩn đoán dùng được: *cùng không gian hoá học + giá trị chung bị
dịch* ⇒ hiệu ứng mẻ thí nghiệm, co giãn lại hấp thụ được (đúng trường hợp clearance);
*khác không gian hoá học + giá trị chung nhất quán* ⇒ khác miền, co giãn lại sẽ phá tín
hiệu (đúng trường hợp half-life).

## 6. Hai phát hiện khác, dùng được ngay

### 6a. AZ clearance cao hơn chuẩn vàng 0,68–0,855 log

Parrondo-Pizarro: TDC có phân bố clearance **cao hơn đáng kể** so với 6 nguồn khác
(KS p < 0,0001), chênh trung bình **0,680–0,855 log đơn vị (4,8–7,2 mL/phút/kg)**. Và
đây **không** phải chuyện miền áp dụng — không gian hoá học chồng lấn mạnh, tức **TDC báo
clearance cao hơn cho chính những phân tử nó chia sẻ với chuẩn vàng**.

Giải thích của họ: TDC là CLint in vitro, bỏ qua gắn protein huyết tương nên **ước lượng
vượt** clearance toàn phần. **Trúng ngay nhánh in vitro của cầu IVIVE của mình** — phải
trích dẫn ở Mục 3.

### 6b. `Half_Life_Obach` được xác nhận là lành

Obach là một trong hai bộ **chuẩn tham chiếu** trong thiết kế của họ. Nguồn bị lệch là
**Fan et al. 2024** (bộ đứng sau ADMETlab 3.0), thấp hơn 0,778–1,026 log đơn vị. Các chất
chung giữa Obach và nguồn khác thì **nhất quán**.

→ Giữ `Half_Life_Obach` làm mục tiêu, và trích dẫn đây làm bằng chứng độc lập rằng nó
lành về phân bố.

## 7. PKSmart — con số định lượng cho nhiễm `/F`

Tính trực tiếp từ file họ phát hành:

- **7/40 dòng có VDss > 10 L/kg = 17,5%**; tất cả đều thuộc phần `source=FDA` (7/42 =
  16,7%). Bốn dòng vượt 20 L/kg. Cao nhất: mobocertinib **58,48 L/kg**.
- **Phần Lombardo sạch**: n=9, VDss cao nhất 7,4 L/kg. → nhiễm **chỉ ở phần suy từ nhãn
  FDA**, tức là hiện vật trích xuất, không phải lựa chọn curation.
- **CL cũng nhiễm, và chưa ai nói**: 5/38 giá trị CL (**13,2%**) **vượt lưu lượng máu gan**
  (~21 mL/phút/kg) — bất khả thi với chất không thải qua thận. Cao nhất melphalan
  flufenamide **192,22**.
- Nhiễm **lan sang** `External_test_315.csv`: đúng bảy giá trị VDss đó xuất hiện lại
  (7/51 = 13,7%).

## Việc tiếp theo

- [ ] Chờ SI của Schuh → giải quyết mâu thuẫn 4,95%.
- [ ] Sửa Mục 5.3 và 5.5 của bài theo phát biểu ở §4.
- [ ] Thêm trích dẫn Parrondo-Pizarro vào Mục 3 (§6a) và Mục 2 (§6b).
- [ ] Khử trùng bộ external của mình với `FDA_v2.csv` (42 tên) và `External_test_315.csv`
      (chỉ có `smiles_r`, phải qua InChIKey, và chuẩn hoá proton trước khi khớp).
