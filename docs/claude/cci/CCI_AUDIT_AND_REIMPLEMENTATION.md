# CCI Audit & Re-implementation Plan

> 2026-02-18 | CellChat + MNN 현재 구현 감사, 문제점 진단, 재구현 완료

## 0. 재구현 결과 요약 (Re-implementation Result)

### CellChat v2 재구현: ✅ 성공
- **전략**: Seurat에서 condition별로 직접 `createCellChat()` 실행 (per-sample merge 우회)
- **핵심 해결**: `@data` slot이 정상적으로 채워져 모든 comparison 함수 작동
- **새로 생성된 플롯 (이전에 MISSING)**:
  - `07_rankNet_stacked.png` / `08_rankNet_unstacked.png` — **Information flow 비교** ✅
  - `09_bubble_comparison.png` — **L-R pair bubble plot** ✅
  - `03/04_diff_circle_{count,weight}.png` — 표준 CellChat diffInteraction ✅
  - `05/06_diff_heatmap_{count,weight}.png` — 표준 CellChat heatmap ✅
  - Per-pathway chord + contribution plots (36 PNG per layer) ✅
- **추가 해결 (fix_cellchat_heatmaps.R)**:
  - `05/06_diff_heatmap_{count,weight}.png` — ComplexHeatmap draw 이슈 해결 ✅
  - `10_outgoing_*.png` / `11_incoming_*.png` — 개별 condition별 signaling role heatmap ✅
  - `12_signaling_scatter.png` + per-celltype scatter — lift→re-merge로 centrality 보존 ✅
  - `13_signaling_role_scatter_*.png` — outgoing vs incoming per cell type ✅
- **L1**: 16 main + 36 pathway + 7 scatter = **59 PNG**
- **L2**: 18 main + 36 pathway + 5 scatter = **59 PNG**
- **Output**: `cci/plots/cellchat_L1_cohort_v2/`, `cellchat_L2_g3_v2/`
- **Scripts**: `scripts/reimplement_cellchat_v2.R`, `scripts/fix_cellchat_v2_final.R`, `scripts/fix_cellchat_heatmaps.R`

### MNN Contrast Bug Fixed & Re-run: ✅ COMPLETE (2026-02-19)
- **이전 오진**: "lr_target_prior_cor sample ≥5 불충족" → 실제로는 contrast string format mismatch
- **진짜 원인**: `limma::makeContrasts(Stroke-HC)` → `"Stroke - HC"` (spaces) vs `contrast_tbl` = `"Stroke-HC"` (no spaces) → `inner_join` 0 rows
- **수정**: `contrasts_oi`를 apostrophe wrapping (`"'Stroke-HC'"`) → `makeContrasts`가 string literal로 처리
- **v3 결과**: L1 anno2: 8,244 / L1 anno1: 148,958 / L2 anno2: 7,600 / L2 anno1: 91,111 interactions
- Output: `cci/mnn/{L1_cohort_anno2_v3, L1_cohort_anno1_v3, L2_g3_anno2_v3, L2_g3_anno1_v3}/`
- 상세: `docs/claude/cci/MNN_FAILURE_ANALYSIS.md`

### L1 rankNet 해석 (새 발견)

**HC-dominant pathways** (HC > Stroke): TGFb, VISFATIN, PARs, TNF, CSF, IL2, IL16
**Stroke-dominant pathways** (Stroke > HC): RESISTIN, BAFF, ANNEXIN, FLT3, MIF, CypA
**균등 pathways**: BAG, GALECTIN

→ SIIS는 전체 CCI 감소이지만, 특정 pro-inflammatory pathways (RESISTIN, MIF, CypA)는 오히려 Stroke에서 증가 — **selective immunosuppression + inflammatory activation** pattern

### L2 rankNet 해석 (새 발견)

**Bad-dominant pathways**: IL1, RESISTIN, FLT3, BAFF, GALECTIN, MIF, CypA, ANNEXIN, IL16
**Good-dominant/balanced**: TNF, BAG

→ g3 Bad에서 거의 모든 pathway가 증가 — **maladaptive broad immune activation** pattern
→ 특히 IL1 (Bad only), RESISTIN (Bad >> Good) = inflammatory damage markers

---

---

## 1. 감사 결과 요약 (Audit Summary)

### CellChat: 부분적으로 정상

| 항목 | 상태 | 비고 |
|------|------|------|
| Per-sample workflow | ✅ OK | createCellChat → identifyOE → computeCommunProb → aggregateNet |
| Assay/normalization | ✅ OK | RNA assay, subsetData() 사용 |
| DB selection | ✅ OK | "Secreted Signaling" only |
| Per-sample 객체 저장 | ✅ OK | 56 samples (L1), 32 samples (L2) |
| Merged 객체 | ⚠️ | @data 비어있음, @net 구조 불명확 |
| Comparison plots | ❌ | mergeCellChat() 실패 → rankNet, bubble, diffInteraction 등 미생성 |
| Summary statistics | ✅ OK | Per-sample count/weight → boxplot 생성됨 |
| Pathway-level analysis | ✅ OK | 대표 sample에서 native plots 생성 |

**핵심 문제**: Condition-level merged 객체의 `@data` slot이 ncol=0. `mergeCellChat()`이 비교 함수에서 실패. 결과적으로 CellChat의 가장 중요한 비교 기능들(rankNet, netVisual_diffInteraction, information flow)을 사용 불가.

### MNN (MultiNicheNet): 실행은 성공했으나 핵심 결과 없음

| 항목 | 상태 | 비고 |
|------|------|------|
| SCE conversion | ✅ OK | RNA assay → SCE |
| Contrast specification | ✅ FIXED | 원래 "Stroke-HC" → apostrophe wrapping "'Stroke-HC'" |
| Cell type filtering | ✅ OK | 18-19/21 cell types retained (pDC, Platelet.PLA excluded) |
| DE analysis | ✅ OK | 83K genes (L1 anno1), 7K (L2 anno2) |
| Ligand activity | ✅ OK | 1.2M (L1 anno1), 106K (L2 anno2) pairs |
| **group_prioritization_tbl** | ✅ FIXED | L1 anno1: 148,958 / L2 anno1: 91,111 interactions |

**이전 문제 (해결됨)**: contrast string format mismatch로 `inner_join` 0 rows. `"Stroke-HC"` (contrast_tbl) vs `"Stroke - HC"` (limma output). Apostrophe wrapping으로 해결. 상세: `MNN_FAILURE_ANALYSIS.md`.

---

## 2. 근본 원인 분석 (Root Cause)

### CellChat

**문제**: `run_cellchat_cli.R`의 merge 단계에서 `mergeCellChat()`이 여러 sample-level 객체를 condition-level로 합칠 때, 각 sample의 cell barcode가 다른 sample과 충돌 → `cell.prefix=TRUE`로 해결하지만, 합쳐진 후 `@data` slot의 column count가 0이 되는 현상.

**근본 원인 가설**:
1. CellChat v2의 merge 로직이 `@data.signaling`만 보존하고 `@data` (raw expression)은 drop
2. 또는 merge 시 barcode prefix 적용이 `@data`와 `@meta`에 일관되지 않음

**영향**: 비교 분석의 핵심인 `rankNet()`, `netVisual_diffInteraction()`, `netVisual_bubble()` 등을 사용할 수 없음. 현재 생성된 비교 플롯은 **custom implementation** (직접 작성한 circle/heatmap/boxplot)이며 CellChat 공식 비교 기능이 아님.

### MNN — RESOLVED ✅ (2026-02-19)

**문제 (해결됨)**: `group_prioritization_tbl = 0 rows` (모든 6회 실행).

**진짜 근본 원인**: Contrast string format mismatch.
- `contrast_tbl$contrast = "Stroke-HC"` (from CLI)
- `celltype_de$contrast = "Stroke - HC"` (from limma::makeContrasts deparse)
- `inner_join` → 0 rows (silent failure, no warning)

**이전 오진**: "lr_target_prior_cor sample ≥5 불충족" — 실제로는 sample 충분 (L1: 17/20 cell types ≥5). 이 메시지는 `group_prioritization_tbl`이 이미 0 rows일 때 `receivers_oi = character(0)`이 되어 출력됨.

**수정**: `contrasts_oi`를 apostrophe wrapping (`"'Stroke-HC'"`) → `makeContrasts`가 string literal로 처리 → column name 보존.

**v3 결과**: L1 anno1 148,958 / L1 anno2 8,244 / L2 anno1 91,111 / L2 anno2 7,600 interactions.

---

## 3. 현재 결과의 신뢰도 (Reliability Assessment)

### 신뢰할 수 있는 결과

1. **CellChat summary statistics** (per-sample interaction count/weight mean): ✅
   - Per-sample CellChat 객체는 정상 → 각 sample의 총 interaction count/weight 합산은 정확
   - L1: HC 114.75 vs Stroke 72.39 (−37%) — 신뢰 가능
   - L2: Good 69.18 vs Bad 76.24 (+10%) — 신뢰 가능

2. **CellChat pathway-level analysis** (대표 sample에서): ⚠️ 제한적
   - 대표 sample 1개의 결과 → population-level inference 불가
   - 어떤 pathway가 존재하는지는 확인 가능, 양적 비교는 부적절

3. **MNN DE genes** (per cell type): ✅
   - multinichenetr의 DE 단계는 muscat-edgeR pseudobulk 사용
   - cDC1 594 sig genes, CD4_S100A8 142 sig genes 등 — 방향성 신뢰 가능
   - 단, 이것은 standard DEG와 동일 (CCI-specific insight 아님)

4. **MNN ligand activities**: ⚠️ 제한적
   - Activity score 자체는 NicheNet prior에 기반한 prediction
   - Experimental validation 없이는 ranking 정도만 참고 가능

### 신뢰할 수 없는 결과

1. **CellChat condition 비교 (circle diff, heatmap diff)**: ❌
   - Custom implementation으로 생성되어 CellChat의 statistical framework 미사용
   - 단순 count/weight 차이를 시각화한 것이지 통계적 비교가 아님

2. **MNN prioritized interactions**: ✅ FIXED (v3)
   - L1 anno1: 148,958 / L2 anno1: 91,111 interactions — 시각화 및 해석 필요

3. **MNN circos/heatmap** (from ligand_activities_target_de_tbl): ⚠️
   - 이 테이블에는 sender 정보가 없음 → ligand→receiver만 표시
   - Prioritization score 없이 raw activity로 ranking → 과해석 위험

---

## 4. 재구현 계획 (Re-implementation Plan)

### Option A: CellChat 재구현 (Recommended)

**목표**: Standard CellChat comparison workflow를 정상 작동시켜 `rankNet()`, `netVisual_diffInteraction()`, `netVisual_bubble()` 등 공식 비교 함수 사용.

**전략**: Per-sample 객체를 직접 로드하여 condition별로 올바르게 merge.

```r
# === CellChat Re-implementation (sound version) ===
library(CellChat)
library(qs)

BASE <- "/data/user3/sobj/stroke_hc_v8_2"

# ---- Step 1: Load per-sample objects ----
sample_dirs <- list.dirs(
  file.path(BASE, "cci/cellchat/L1_cohort_anno2/samples"),
  recursive = FALSE, full.names = TRUE
)
cc_list <- lapply(sample_dirs, function(d) {
  cc <- qread(file.path(d, "cellchat.qs"))
  cc <- updateCellChat(cc)  # v2 compatibility
  return(cc)
})
names(cc_list) <- basename(sample_dirs)

# ---- Step 2: Map samples to conditions ----
# Load Seurat metadata to get patient → condition mapping
sobj <- qread(file.path(BASE, "5_1_hc_is.qs"))
mapping <- unique(sobj@meta.data[, c("patient_name", "cohort")])
mapping$dir_name <- make.names(mapping$patient_name)

# Split sample list by condition
hc_names <- mapping$dir_name[mapping$cohort == "HC"]
is_names <- mapping$dir_name[mapping$cohort != "HC"]

cc_hc <- cc_list[names(cc_list) %in% hc_names]
cc_is <- cc_list[names(cc_list) %in% is_names]
rm(sobj); gc()

# ---- Step 3: Merge within condition (keeping per-sample structure) ----
cc_hc_merged <- mergeCellChat(cc_hc, add.names = names(cc_hc))
cc_is_merged <- mergeCellChat(cc_is, add.names = names(cc_is))

# ---- Step 4: Create comparison list ----
cc_compare <- list(HC = cc_hc_merged, Stroke = cc_is_merged)
cellchat <- mergeCellChat(cc_compare, add.names = c("HC", "Stroke"))

# ---- Step 5: Standard comparison functions ----
# Information flow
rankNet(cellchat, mode = "comparison")

# Differential interaction
netVisual_diffInteraction(cellchat, weight.scale = TRUE)

# Bubble plot comparison
netVisual_bubble(cellchat, comparison = c(1, 2), remove.isolate = FALSE)

# Heatmap
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", signaling = NULL)
```

**주의사항**:
- Per-sample merge 후 condition merge를 2단계로 수행
- `updateCellChat()` 필수 (CellChat v2 호환성)
- Cell type가 모든 sample에 존재하지 않을 수 있음 → `filterCommunication()`에서 min.cells로 조절

**대안 (simpler)**: Per-sample objects를 condition별로 바로 merge하지 않고, **liftCellChat()**으로 cell type space를 통일한 후 merge:
```r
# Union of all cell types across all samples
all_ct <- sort(unique(unlist(lapply(cc_list, function(x) levels(x@idents)))))
cc_list <- lapply(cc_list, function(x) liftCellChat(x, group.new = all_ct))
```

### Option B: MNN 재구현 (with relaxed thresholds)

**목표**: `group_prioritization_tbl`에 실제 결과 확보.

**전략 1: 더 관대한 parameters**

```r
Rscript "$MNN_SCRIPT" \
  -i "$BASE/5_1_hc_is.qs" \
  -g cohort \
  -s patient_name \
  -c anno1 \
  -f "Stroke-HC" \
  --min_cells 5 \            # 10 → 5
  --logfc_threshold 0.05 \   # 0.10 → 0.05
  --fraction_cutoff 0.01 \   # 0.05 → 0.01
  --p_val_threshold 0.10 \   # 0.05 → 0.10
  --cores 16 \
  -o "${BASE}/cci/mnn/L1_cohort_anno1_permissive"
```

**전략 2: anno2 level에서 재시도**

anno1 (19 types)은 combination이 너무 많아 sparse. anno2 (8 compartments)로 돌아가되, thresholds를 더 permissive하게:
```r
--min_cells 5 --logfc_threshold 0.05 --fraction_cutoff 0.01 --p_val_threshold 0.10
```

**전략 3: Custom prioritization**

MNN의 intermediate outputs (`de_genes_df`, `ligand_activities`)를 활용하여 직접 prioritization:

```r
# Load MNN results
mnn <- qread("multinichenet_results.qs")

# Extract DE genes as potential receptors/ligands
de_genes <- mnn$ligand_activities_targets_DEgenes$de_genes_df
la <- mnn$ligand_activities_targets_DEgenes$ligand_activities

# Load LR network
lr_network <- readRDS("/data/user3/git_repo/human/lr_network_human_21122021.rds")

# Manual prioritization: find DE genes that are also ligands/receptors
de_ligands <- de_genes %>%
  filter(p_adj < 0.05) %>%
  inner_join(lr_network, by = c("gene" = "from"))  # gene is ligand

de_receptors <- de_genes %>%
  filter(p_adj < 0.05) %>%
  inner_join(lr_network, by = c("gene" = "to"))    # gene is receptor

# Build sender-receiver-ligand-receptor table manually
# ...
```

### Option C: Alternative CCI Method — LIANA

**목표**: CellChat/MNN이 문제가 많으면, LIANA (LIgand-receptor ANAlysis)로 대체.

LIANA는 여러 CCI method (CellChat, CellPhoneDB, NATMI, Connectome, SingleCellSignalR)의 consensus를 제공하며, condition comparison이 내장되어 있음.

```r
# R implementation
library(liana)

# Basic run (multiple methods, consensus)
liana_res <- liana_wrap(sce, method = c("natmi", "connectome", "sca", "cellchat"))
liana_agg <- liana_aggregate(liana_res)

# Condition-specific (LIANA+)
# Differential CCI between conditions
liana_diff <- liana_bycluster(
  sce = sce,
  group_col = "cohort",
  sample_col = "patient_name",
  celltype_col = "anno1"
)
```

**장점**:
- Multi-method consensus → single method bias 방지
- Built-in condition comparison
- Better maintained than multinichenetr
- 논문에서 "multi-method CCI consensus" approach로 제시 가능

---

## 5. 권장 실행 순서 (Recommended Action)

```
Priority 1 (즉시):
├── Option A: CellChat 재구현 — per-sample → liftCellChat → condition merge → 비교 함수
│   - 기존 per-sample 객체 활용 가능 (재실행 불필요)
│   - 예상 소요: 1-2시간 (스크립트 작성 + 실행)
│
Priority 2 (이후):
├── Option B-1: MNN permissive re-run (anno1, relaxed thresholds)
│   - 예상 소요: 30-60분 실행
│
Priority 3 (대안):
├── Option C: LIANA 설치 + 실행
│   - CellChat/MNN 모두 문제 시 대안
│   - 예상 소요: 설치 1시간 + 실행 2-3시간
```

---

## 6. MNN Wrapper 수정 필요 사항

파일: `myR/R/cci_multinichenet_wrapper.R`

1. **Line 53**: 중복 `p_val_thresh = 0.05` 제거
2. **Line 28**: 중복 `@importFrom tibble tibble as_tibble` 제거
3. **Line 30-32**: 하드코딩된 source() 경로 제거 또는 조건부 로드
4. **결과 저장**: `sender_receiver_de`, `sender_receiver_info` 등 intermediate tables도 저장하도록 추가

---

## 7. docs/claude/ 에서 참조 가능한 스크립트

| Script | Location | Purpose |
|--------|----------|---------|
| CCI 전체 실행 | `results/scripts/run_cci_all.sh` | CellChat + MNN dual-layer |
| MNN anno1 실행 | `results/scripts/run_mnn_anno1.sh` | MNN re-run at anno1 |
| CellChat native plots | `results/scripts/plot_cci_comprehensive.R` | Pathway circles, bubbles |
| MNN anno1 plots | `results/scripts/plot_mnn_anno1.R` | DE volcanos, circos, heatmaps |
| CellChat CLI | `Git_Repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R` | Main CellChat pipeline |
| MNN CLI | `Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R` | Main MNN pipeline |
| MNN wrapper | `myR/R/cci_multinichenet_wrapper.R` | myR package wrapper |
