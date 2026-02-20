# MultiNicheNet (MNN) Failure Analysis: `group_prioritization_tbl` = 0 Rows

> 2026-02-19 | 전체 6회 실행 모두 0 rows인 원인 분석 및 해결 방안

---

## 1. 문제 (Problem)

MultiNicheNet 분석의 최종 핵심 출력인 `group_prioritization_tbl`이 **모든 실행에서 0 rows**:

| Run | Layer | Cell Types | Thresholds | DE genes | Ligand Activities | **Prioritized** |
|-----|-------|-----------|------------|----------|-------------------|-----------------|
| L1_cohort_anno2 | L1 (HC vs Stroke) | anno2 (6) | default | 5,511 | 231,116 | **0** |
| L1_cohort_anno1 | L1 (HC vs Stroke) | anno1 (19) | default | 78,619 | 1,237,237 | **0** |
| L1_cohort_anno2_permissive | L1 (HC vs Stroke) | anno2 (6) | relaxed | 16,729 | 378,508 | **0** |
| L2_g3_anno2 | L2 (g3=1 vs g3=2) | anno2 (6) | default | 6,186 | 106,053 | **0** |
| L2_g3_anno1 | L2 (g3=1 vs g3=2) | anno1 (18) | default | 70,095 | 525,159 | **0** |
| L2_g3_anno2_permissive | L2 (g3=1 vs g3=2) | anno2 (6) | relaxed | 19,324 | 250,748 | **0** |

주목: DE analysis, ligand activity 계산 등 **중간 단계는 모두 정상**. 문제는 최종 prioritization 단계에서만 발생.

---

## 2. 근본 원인 (Root Cause): Contrast String Format Mismatch

### 진단 결과

`group_prioritization_tbl`이 0 rows인 **유일한 원인**은 `contrast` 컬럼의 문자열 불일치:

```
contrast_tbl$contrast  = "Stroke-HC"      ← CLI에서 전달된 원본 (공백 없음)
celltype_de$contrast   = "Stroke - HC"    ← muscat/limma가 생성한 값 (공백 있음)
```

`generate_prioritization_tables()` 함수 내부에서 가장 첫 번째 `inner_join`에서 실패:

```r
# generate_prioritization_tables() 내부:
group_prioritization_tbl = contrast_tbl %>%
  inner_join(sender_receiver_de) %>%    # ← 여기서 0 rows (contrast 불일치!)
  inner_join(ligand_activities) %>%     # ← 입력이 0이므로 0
  inner_join(avg_df_group) %>%          # ← 입력이 0이므로 0
  ...
```

### 원인 체인 (Causation Chain)

```
1. User CLI:  -f "Stroke-HC" (공백 없음)
       ↓
2. Wrapper (cci_multinichenet_wrapper.R, line 141-143):
   contrast_tbl <- tibble(contrast = contrasts_oi, group = groups)
   → contrast_tbl$contrast = "Stroke-HC"
       ↓
3. multi_nichenet_analysis() → perform_muscat_de_analysis():
   limma::makeContrasts(Stroke-HC, levels=design)
   → contrast matrix column name = "Stroke - HC" (limma adds spaces around operators)
       ↓
4. muscat::resDS() → de_output_tidy$contrast = "Stroke - HC"
       ↓
5. combine_sender_receiver_de() inherits contrast from celltype_de
   → sender_receiver_de$contrast = "Stroke - HC"
       ↓
6. generate_prioritization_tables():
   contrast_tbl (contrast="Stroke-HC") %>%
     inner_join(sender_receiver_de (contrast="Stroke - HC"))
   → 0 rows!!!
```

### 증명 (Verification)

R에서 직접 확인:

```r
# L1_cohort_anno2_permissive 결과에서:
contrast_tbl_wrong  <- tibble(contrast = "Stroke-HC", group = "Stroke")
contrast_tbl_correct <- tibble(contrast = "Stroke - HC", group = "Stroke")

contrast_tbl_wrong  %>% inner_join(sender_receiver_de)  # → 0 rows
contrast_tbl_correct %>% inner_join(sender_receiver_de)  # → 2,524 rows
```

L2도 동일:
```
contrast_tbl$contrast = "X2-X1"   vs   celltype_de$contrast = "X2 - X1"
```

### 이 버그는 어디서 발생하는가?

**multinichenetr 패키지 자체의 설계 문제**. `multi_nichenet_analysis()` 함수는:
1. `contrast_tbl`을 사용자로부터 받아 그대로 유지
2. `perform_muscat_de_analysis()`에서 `limma::makeContrasts()`를 호출하여 DE를 수행
3. `limma::makeContrasts()`는 `"A-B"` 형식을 `"A - B"` (공백 포함) 형식으로 변환
4. `muscat::resDS()`는 이 변환된 이름을 `contrast` 컬럼에 사용
5. `generate_prioritization_tables()`에서 두 형식을 `inner_join`으로 매칭 시도

다시 말해, **multinichenetr는 사용자가 `contrast_tbl$contrast`를 limma 형식 (`"A - B"`)으로 제공한다고 가정**하지만, 이를 문서화하지 않았고 validation도 하지 않음. 또한 CLI wrapper에서 `-f "Stroke-HC"` 형식을 그대로 `contrast_tbl`에 넣으면서 mismatch 발생.

---

## 3. 오진 (Misdiagnosis): "lr_target_prior_cor" 문제

### 로그에서 본 메시지

모든 6회 실행에서 다음 메시지가 출력됨:

```
"For no celltypes, sufficient samples (>= 5) were available for a correlation analysis.
 lr_target_prior_cor, the output of this function, will be NULL."
```

### 왜 이것이 오해를 야기했는가

이 메시지로 인해 다음과 같이 진단됨:
> "근본 원인: lr_target_prior_cor step에서 cell type × condition당 sample ≥ 5 필요 → 불충족"

그러나 **이것은 진짜 원인이 아니다**. `lr_target_prior_cor_inference()` 함수는:

```r
lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
  ...)
```

`group_prioritization_tbl`이 이미 0 rows이므로, `receivers_oi = character(0)`. `lapply(character(0), ...)` 는 빈 list를 반환하고, 함수는 "sufficient samples" 메시지와 함께 빈 tibble을 반환한다. **실제 sample 수와 무관하게** 이 메시지가 출력된다.

### 실제 sample 수 (충분함)

**L1 (HC vs Stroke, anno2, min_cells=5):**

| Cell Type | HC Samples | Stroke Samples | ≥5 Both? |
|-----------|-----------|----------------|----------|
| Bc | 20 | 36 | Yes |
| Tc | 20 | 36 | Yes |
| Mono | 20 | 36 | Yes |
| NKc | 20 | 36 | Yes |
| DC | 20 | 33 | Yes |
| Platelet/PLA | 20 | 4 | **No** |

**5/6 cell types가 ≥5 samples in both conditions**. Sample 수는 충분했다.

**L1 (HC vs Stroke, anno1, min_cells=10):**

| Cell Type | HC Samples | Stroke Samples | ≥5 Both? |
|-----------|-----------|----------------|----------|
| B_cell | 20 | 36 | Yes |
| CD4+ T_Naive/Memory | 20 | 36 | Yes |
| CD8+ T_Cytotoxic | 20 | 36 | Yes |
| Inflammatory Monocyte | 20 | 36 | Yes |
| NK_cell | 20 | 36 | Yes |
| MAIT | 20 | 35 | Yes |
| Plasma_cell | 20 | 35 | Yes |
| Treg | 20 | 35 | Yes |
| CD16+ Monocyte | 20 | 33 | Yes |
| CD14+ Monocyte | 20 | 32 | Yes |
| CD4_S100A8_CD14 | 20 | 32 | Yes |
| CD8+ Trm | 20 | 32 | Yes |
| ISG+ Myeloid | 20 | 32 | Yes |
| CD4_S100A8_CSF3R | 18 | 23 | Yes |
| cDC2 | 20 | 23 | Yes |
| ISG+ T_cell | 20 | 10 | Yes |
| Platelet/PLA | 18 | 0 | **No** |
| pDC | 17 | 1 | **No** |
| Proliferating | 16 | 7 | Yes |
| cDC1 | 4 | 2 | **No** |

**17/20 cell types가 ≥5 samples in both conditions**. Correlation step은 문제없이 작동해야 함.

**L2 (g3=1 vs g3=2, anno2, min_cells=5):**

| Cell Type | g3=1 Samples | g3=2 Samples | ≥5 Both? |
|-----------|-------------|-------------|----------|
| Bc | 11 | 21 | Yes |
| Tc | 11 | 21 | Yes |
| Mono | 11 | 21 | Yes |
| NKc | 11 | 21 | Yes |
| DC | 10 | 20 | Yes |
| Platelet/PLA | 1 | 3 | **No** |

**5/6 cell types가 ≥5 samples in both conditions**.

**L2 (g3=1 vs g3=2, anno1, min_cells=10):**

15/19 cell types with ≥5 samples in both conditions. (cDC1, pDC, Platelet/PLA, ISG+ T_cell only failed.)

---

## 4. 해결 시도 (Attempted Fixes)

### 시도 1: Relaxed Thresholds (2026-02-18)

**Script**: `/data/user3/sobj/stroke_hc_v8_2/scripts/reimplement_mnn_permissive.sh`

Parameters relaxed:
- `min_cells`: 10 → 5
- `logfc_threshold`: 0.10 → 0.05
- `fraction_cutoff`: 0.05 → 0.01
- `p_val_threshold`: 0.05 → 0.10

**결과**: 4개 runs 모두 0 rows. Threshold 완화는 DE genes 수를 증가시켰지만 (5,511 → 16,729 for L1 anno2), contrast string mismatch는 threshold과 무관하므로 효과 없음.

### 시도 2: anno2 Level (Fewer Cell Types = More Samples Per Type)

**가설**: anno1 (19 types)에서 cell type당 sample이 부족 → anno2 (8 types)로 aggregation하면 sample 증가

**결과**: 0 rows. Sample 수는 실제로 충분했고, 문제는 sample 수가 아님.

### 시도 3: anno1 with Default/Permissive

4가지 조합을 모두 시도:
- L1 × anno1 × default: 0 rows
- L1 × anno2 × permissive: 0 rows
- L2 × anno1 × default: 0 rows
- L2 × anno2 × permissive: 0 rows

### 왜 모든 시도가 실패했는가

**모든 시도가 동일한 CLI wrapper를 사용**하여 `contrast_tbl$contrast`가 항상 공백 없는 형식 (`"Stroke-HC"` 또는 `"X2-X1"`)으로 설정됨. Threshold, cell type resolution, sample 수와 무관하게 contrast string mismatch가 존재하는 한 항상 0 rows.

---

## 5. 왜 이전 진단이 틀렸는가 (Why Previous Diagnosis Was Wrong)

### 이전 진단 (CCI_AUDIT_AND_REIMPLEMENTATION.md)

> "근본 원인: multinichenetr의 correlation 단계는 sender-receiver 쌍별로 최소 5개 sample이 필요.
> anno1 (19 cell types) × 2 conditions에서 sparse cell types (cDC1, ISG+ T_cell 등)는
> condition당 2-3 sample에서만 충분한 cell 수를 가짐."

### 실제 상황

1. **Sample 수는 대부분 충분**: L1 anno2에서 5/6, L1 anno1에서 17/20, L2 anno2에서 5/6 cell types가 ≥5 samples in both conditions
2. **"sufficient samples" 메시지는 오해의 소지**: 이 메시지는 `lr_target_prior_cor_inference()` 함수에서 출력되지만, 실제로는 `group_prioritization_tbl`이 이미 0 rows이므로 `receivers_oi = character(0)`이 되어 correlation 분석 자체가 실행되지 않은 것
3. **진짜 문제는 그 이전 단계**: `generate_prioritization_tables()`의 첫 번째 `inner_join`에서 contrast string 불일치로 0 rows

### 오진의 원인

1. Log 메시지가 `lr_target_prior_cor` 문제를 가리켜서 correlation step에 집중
2. Multinichenetr 코드가 `inner_join`을 silent하게 사용 (0 rows여도 error/warning 없음)
3. Sample 수 부족이 sparse cell types에서 실제로 존재하기는 하므로 (cDC1, pDC 등) 가설이 그럴듯해 보임
4. `group_prioritization_tbl`과 `lr_target_prior_cor` 모두 0 rows여서 correlation 문제처럼 보임

---

## 6. 수정 방법 (Fix) — APPLIED & VERIFIED ✅

### 최종 해결: Apostrophe-wrapping

`myR/R/cci_multinichenet_wrapper.R`에서 `contrasts_oi`를 apostrophe로 감싸기:

```r
# BEFORE (bug): contrasts_oi = "Stroke-HC"
#   → makeContrasts(Stroke-HC, levels=design)
#   → column name = "Stroke - HC" (spaces added by deparse)
#   → contrast_tbl$contrast = "Stroke-HC" ≠ celltype_de$contrast = "Stroke - HC"

# AFTER (fix): contrasts_oi = "'Stroke-HC'"
#   → makeContrasts("Stroke-HC", levels=design)
#   → column name = "Stroke-HC" (preserved as string literal)
#   → contrast_tbl$contrast = "Stroke-HC" == celltype_de$contrast = "Stroke-HC" ✅
```

핵심 원리: `limma::makeContrasts()`는 bare expression `Stroke-HC`를 `"Stroke - HC"`로 deparse하지만,
quoted string `"Stroke-HC"` (apostrophe-wrapped)는 문자열 리터럴로 취급하여 원본을 보존함.

### 수정된 코드 (`cci_multinichenet_wrapper.R`)

```r
# Parse raw contrast names (strip any existing apostrophes)
raw_contrasts <- gsub("'", "", unlist(strsplit(paste(contrasts_oi, collapse = ","), ",")))
raw_contrasts <- trimws(raw_contrasts)
raw_contrasts <- raw_contrasts[raw_contrasts != ""]

# Extract the "group of interest" = first condition in each contrast
groups <- sapply(raw_contrasts, function(x) strsplit(x, "-")[[1]][1])

# Build contrast_tbl: contrast names WITHOUT apostrophes or spaces
contrast_tbl <- tibble::tibble(
    contrast = raw_contrasts,
    group = groups
)

# CRITICAL: wrap in apostrophes for makeContrasts to preserve as string literal
contrasts_oi <- paste0("'", raw_contrasts, "'", collapse = ",")
```

### 검증 결과

L1 anno2 (HC vs Stroke) 테스트:
- **이전**: 0 rows (모든 6회 실행)
- **수정 후**: **8,244 prioritized interactions** ✅
  - 5 senders, 5 receivers, 331 ligands, 337 receptors
  - Top: Mono→Mono S100A8-CD36 (score=0.978)
  - 실행 시간: 7.7분

### 수정 파일

| File | Status |
|------|--------|
| `myR/R/cci_multinichenet_wrapper.R` | ✅ Fixed |
| `/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci_multinichenet_wrapper.R` | ✅ Fixed (same change) |
| CLI script `run_multinichenet.R` | No change needed (wrapper handles apostrophes) |

---

## 7. 현재 활용 방안 (Current Workaround)

`group_prioritization_tbl`은 비어있지만, **중간 결과물은 모두 유효**:

### 활용 가능한 중간 결과

| Output | Content | Rows (L1 anno1) | Usability |
|--------|---------|-----------------|-----------|
| `celltype_de` | Per-cell-type DE genes | 78,619 | Standard muscat-edgeR pseudobulk DE (CCI-specific이 아님) |
| `ligand_activities` | NicheNet ligand activity scores | 1,237,237 | NicheNet prior에 기반한 prediction; ranking 참고용 |
| `ligand_activities_target_de_tbl` | Ligand → target DE 연결 | 378,500 (L1 anno2 permissive) | Sender 정보 없음, ligand→receiver만 |
| `celltype_info` (pb/avg/frq) | Pseudobulk expression | Full | 정상적인 expression data |

### Manual Concordant Join

기존 중간 결과를 활용한 수동 prioritization:

```r
# 1. DE genes 중 ligand/receptor인 것 추출
de_ligands <- celltype_de %>%
  filter(p_val < 0.05, abs(logFC) > 0.1) %>%
  inner_join(lr_network, by = c("gene" = "from")) %>%
  rename(ligand = gene, receptor = to, sender = cluster_id)

de_receptors <- celltype_de %>%
  filter(p_val < 0.05, abs(logFC) > 0.1) %>%
  inner_join(lr_network, by = c("gene" = "to")) %>%
  rename(receptor = gene, ligand = from, receiver = cluster_id)

# 2. L-R pair 구성 (sender DE ligand + receiver DE receptor)
lr_pairs <- de_ligands %>%
  inner_join(de_receptors, by = c("ligand", "receptor"), suffix = c("_sender", "_receiver"))

# 3. Ligand activity score 추가
lr_prioritized <- lr_pairs %>%
  inner_join(ligand_activities, by = c("ligand", "receiver" = "receiver"))
```

---

## 8. 현재 상태 (Current Status) — RESOLVED ✅

### MNN 수정 후 재실행 완료

**Wrapper 코드 수정 + 재실행 결과:**

| Config | Cell Types | Prioritized Interactions | Time |
|--------|-----------|-------------------------|------|
| L1 anno2 (HC vs Stroke) | 5 | **8,244** | 7.7 min |
| L1 anno1 (HC vs Stroke) | 18 | **148,958** | 37.8 min |
| L2 anno2 (g3=1 vs g3=2) | 5 | **7,600** | 4.5 min |
| L2 anno1 (g3=1 vs g3=2) | 17 | **91,111** | ~18 min |
| **Total** | | **255,913** | |

### 보완 분석: CellChat (이미 완료)

CellChat v2 re-implementation 성공적으로 완료. `rankNet()`, `netVisual_diffInteraction()`, `netVisual_bubble()` 등 모든 비교 함수 작동 중.

**Output**: `results/cci/plots/cellchat_{L1_cohort,L2_g3}_v2/` (59 PNG per layer)

### 전략: CellChat (gross, anno2) + MNN (fine, anno1)

CellChat로 compartment-level (anno2) gross interaction 패턴을 보여주고,
MNN으로 fine cell type (anno1) level의 specific L-R interaction을 보여주는 dual approach.

---

## 9. 실행 기록 (Run History)

### 전체 실행 기록

| Date | Run | Script | Duration | Result |
|------|-----|--------|----------|--------|
| 2026-02-12 | L1_cohort_anno2 | run_cci_all.sh | ~9 min | 0 rows |
| 2026-02-12 | L2_g3_anno2 | run_cci_all.sh | ~5 min | 0 rows |
| 2026-02-16 | L1_cohort_anno1 | run_mnn_anno1.sh | ~35 min | 0 rows |
| 2026-02-16 | L2_g3_anno1 | run_mnn_anno1.sh | ~15 min | 0 rows |
| 2026-02-18 | L1_cohort_anno2_permissive | reimplement_mnn_permissive.sh | ~9 min | 0 rows |
| 2026-02-18 | L2_g3_anno2_permissive | reimplement_mnn_permissive.sh | ~5 min | 0 rows |
| **2026-02-19** | **L1_cohort_anno2_v3** | **rerun_mnn_fixed.sh** | **7.7 min** | **8,244 rows ✅** |
| **2026-02-19** | **L1_cohort_anno1_v3** | **rerun_mnn_all_v3.sh** | **37.8 min** | **148,958 rows ✅** |
| **2026-02-19** | **L2_g3_anno2_v3** | **rerun_mnn_all_v3.sh** | **4.5 min** | **7,600 rows ✅** |
| **2026-02-19** | **L2_g3_anno1_v3** | **rerun_mnn_all_v3.sh** | **~18 min** | **91,111 rows ✅** |

### 결과 파일 위치

```
/data/user3/sobj/stroke_hc_v8_2/cci/mnn/
├── L1_cohort_anno1/
│   ├── multinichenet_results.qs        # 78K DE genes, 1.2M ligand activities, 0 prioritized
│   ├── multinichenet_params.txt
│   └── run_multinichenet.log
├── L1_cohort_anno2/
│   ├── multinichenet_results.qs        # 5.5K DE genes, 231K ligand activities, 0 prioritized
│   └── ...
├── L1_cohort_anno2_permissive/
│   ├── multinichenet_results.qs        # 16.7K DE genes, 379K ligand activities, 0 prioritized
│   └── ...
├── L2_g3_anno1/
│   ├── multinichenet_results.qs        # 70K DE genes, 525K ligand activities, 0 prioritized
│   └── ...
├── L2_g3_anno2/
│   ├── multinichenet_results.qs        # 6.2K DE genes, 106K ligand activities, 0 prioritized
│   └── ...
├── L2_g3_anno2_permissive/
│   ├── multinichenet_results.qs        # 19.3K DE genes, 251K ligand activities, 0 prioritized
│   └── ...
├── L1_cohort_anno2_v3/                 # ✅ FIXED
│   └── multinichenet_results.qs        # 8,244 prioritized interactions
├── L1_cohort_anno1_v3/                 # ✅ FIXED
│   └── multinichenet_results.qs        # 148,958 prioritized interactions
├── L2_g3_anno2_v3/                     # ✅ FIXED
│   └── multinichenet_results.qs        # 7,600 prioritized interactions
└── L2_g3_anno1_v3/                     # ✅ FIXED
    └── multinichenet_results.qs        # 91,111 prioritized interactions
```

### 로그 파일

```
/data/user3/sobj/stroke_hc_v8_2/logs/
├── mnn_L1_cohort_anno2.log
├── mnn_L2_g3_anno2.log
├── mnn_anno1_master.log
├── mnn_L1_anno2_permissive.log
├── mnn_L2_anno2_permissive.log
└── mnn_plots.log
```

---

## 10. 관련 코드 위치 (Code References)

### Bug Location

| File | Line | Issue |
|------|------|-------|
| `myR/R/cci_multinichenet_wrapper.R` | 141-143 | `contrast_tbl` construction without limma format |
| multinichenetr::generate_prioritization_tables | First `inner_join` | Silent 0-row join on `contrast` mismatch |
| multinichenetr::perform_muscat_de_analysis | 267-268 | `limma::makeContrasts()` adds spaces |
| multinichenetr::lr_target_prior_cor_inference | 140-141 | Misleading "insufficient samples" message |

### Key Function Call Chain

```
run_multinichenet.R (CLI)
  → run_multinichenet_analysis() [wrapper]
    → contrast_tbl = tibble(contrast="Stroke-HC", ...)
    → multi_nichenet_analysis() [multinichenetr]
      → perform_muscat_de_analysis()
        → limma::makeContrasts(Stroke-HC, levels=design)
          → column name = "Stroke - HC"  ← SPACES ADDED
        → muscat::resDS()
          → de_output_tidy$contrast = "Stroke - HC"
      → celltype_de$contrast = "Stroke - HC"
      → combine_sender_receiver_de()
        → sender_receiver_de$contrast = "Stroke - HC"
      → generate_prioritization_tables(contrast_tbl, sender_receiver_de, ...)
        → contrast_tbl$contrast ("Stroke-HC") ≠ sender_receiver_de$contrast ("Stroke - HC")
        → inner_join → 0 rows
      → lr_target_prior_cor_inference(receivers_oi = character(0), ...)
        → "For no celltypes, sufficient samples >= 5..."  ← MISLEADING
```

### multinichenetr Package Info

```
Package: multinichenetr
Location: /home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu/multinichenetr
State: compiled (R/multinichenetr.rdb), source not editable
```

---

## 11. 교훈 (Lessons Learned)

1. **Silent `inner_join` failures**: `dplyr::inner_join`은 0 rows를 반환해도 error/warning이 없음. 대규모 pipeline에서 중간 단계의 row count를 반드시 확인해야 함.

2. **Misleading error messages**: `lr_target_prior_cor_inference`의 "insufficient samples" 메시지가 실제 원인을 숨김. 이 함수는 `group_prioritization_tbl`이 0 rows일 때도 동일한 메시지를 출력.

3. **Format consistency**: `limma::makeContrasts()`는 수식의 연산자 주위에 공백을 추가함 (`"A-B"` → `"A - B"`). 이 변환은 muscat의 DE output에 전파되며, contrast_tbl을 수동으로 구성할 때 이를 반영해야 함.

4. **Debug priority**: Intermediate table의 row count를 확인하는 것이 parameter tuning보다 선행되어야 함. 6번의 재실행이 모두 동일한 버그 때문에 실패했음.

---

## 부록 A: multinichenetr `generate_prioritization_tables` 내부 Join Chain

```r
generate_prioritization_tables <- function(
  sender_receiver_info,     # avg_df_group, frq_df_group, pb_df_group (with sender/receiver/ligand/receptor)
  sender_receiver_de,       # contrast/sender/receiver/ligand/receptor/lfc_*/p_val_*
  ligand_activities_targets_DEgenes,  # contrast/receiver/ligand/activity/activity_scaled
  contrast_tbl,             # contrast/group  ← BUG: "Stroke-HC" instead of "Stroke - HC"
  sender_receiver_tbl,      # sender/receiver
  grouping_tbl,             # sample/group
  ...
) {
  # Step 1: Prioritize receptor DE
  receiver_receptor_prioritization = sender_receiver_de %>% ...  # Uses sender_receiver_de directly

  # Step 2: Prioritize ligand activity (up & down)
  receiver_ligand_activity_prioritization_up = ligand_activities %>% filter(direction_regulation == "up") %>% ...

  # Step 3: Prioritize sender ligand DE
  sender_ligand_prioritization = sender_receiver_de %>% ...

  # Step 4: Expression specificity from sender_receiver_info$avg_df_group
  ligand_celltype_specificity = sender_receiver_info$avg_df_group %>%
    inner_join(sender_receiver_tbl) %>%  # join on sender
    inner_join(contrast_tbl) %>%         # join on group  ← THIS WORKS (group matches)
    ...

  # Step 5: L-R fraction
  ligand_receptor_expressed = sender_receiver_info$frq_df %>%
    inner_join(grouping_tbl) %>% ...

  # CRITICAL Step 6: Build final table
  group_prioritization_tbl = contrast_tbl %>%     # contrast="Stroke-HC", group="Stroke"
    inner_join(sender_receiver_de) %>%            # ← FAIL: contrast mismatch!
    inner_join(ligand_activities) %>%
    inner_join(avg_df_group) %>%
    inner_join(frq_df_group) %>%
    inner_join(sender_ligand_prioritization) %>%
    inner_join(receiver_receptor_prioritization) %>%
    inner_join(receiver_ligand_activity_prioritization_up) %>%
    inner_join(receiver_ligand_activity_prioritization_down) %>%
    inner_join(ligand_celltype_specificity) %>%
    inner_join(receptor_celltype_specificity) %>%
    inner_join(ligand_receptor_expressed) %>%
    mutate(prioritization_score = ...)
}
```

---

## 부록 B: Wrapper 수정 (Bug Fix) — APPLIED ✅

**파일**: `myR/R/cci_multinichenet_wrapper.R` (양쪽 copy 모두 수정)

**이전 (BUG)**: contrasts_oi를 그대로 contrast_tbl에 넣고, makeContrasts에 bare expression으로 전달
```r
groups <- sapply(contrasts_oi, function(x) strsplit(x, "-")[[1]][1])
contrast_tbl <- tibble::tibble(contrast = contrasts_oi, group = groups)  # "Stroke-HC"
# → makeContrasts(Stroke-HC, levels=design) → deparse → "Stroke - HC" → MISMATCH
```

**시도 1 (실패)**: gsub로 spaces 추가 → `contrasts_oi = "Stroke - HC"`
→ `multi_nichenet_analysis` 내부 validation에서 condition 추출 실패

**시도 2 (실패)**: contrasts_oi에만 spaces 추가, contrast_tbl은 원본 유지
→ validation이 `contrasts_oi %in% contrast_tbl$contrast` 확인 → 불일치

**최종 수정 (SUCCESS ✅)**: Apostrophe-wrapping
```r
raw_contrasts <- gsub("'", "", unlist(strsplit(paste(contrasts_oi, collapse = ","), ",")))
raw_contrasts <- trimws(raw_contrasts); raw_contrasts <- raw_contrasts[raw_contrasts != ""]
groups <- sapply(raw_contrasts, function(x) strsplit(x, "-")[[1]][1])
contrast_tbl <- tibble::tibble(contrast = raw_contrasts, group = groups)  # "Stroke-HC"
contrasts_oi <- paste0("'", raw_contrasts, "'", collapse = ",")  # "'Stroke-HC'"
# → makeContrasts("Stroke-HC", levels=design) → string literal → "Stroke-HC" preserved ✅
```
