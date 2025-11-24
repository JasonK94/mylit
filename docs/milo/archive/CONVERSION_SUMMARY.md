# TNBC Myeloid h5ad → Seurat 변환 & Milo 파이프라인 정리

본 문서는 `/data/ARPAH/250305_TNBC_26/.../myeloid_v1_TNBC.h5ad`를 Seurat 객체(.qs)로 변환하고 Milo 분석까지 이어가는 과정에서 구현/검증한 내용을 정리한 것이다. 변환 모듈 확장, 변환 스크립트, 산출물 구조, 워크트리 구성 등을 한 번에 파악할 수 있도록 구성했다.

---

## 1. 변환 모듈(`myR/R/h5ad_to_seurat.R`) 핵심 변경 사항

| 기능 | 설명 |
| --- | --- |
| `required_meta_cols` | 필요한 obs 컬럼만 선택하여 저장 (기본: 전체) |
| `required_obsm_keys` | 필요한 obsm만 Seurat reduction으로 변환 (기본: 전체) |
| `store_var` | `adata$var`를 가능한 경우 assay meta.features로, 실패 시 `@misc$anndata_var`로 저장 |
| `store_uns` | `adata$uns`를 `@misc$anndata_uns`에 저장 |
| `store_obsp` | `adata$obsp`가 존재하면 `@misc$anndata_obsp`에 sparse matrix로 저장 |

> TNBC Myeloid h5ad에는 `obsp`가 없기 때문에 저장 시도만 하고 넘어간다.

---

## 2. 변환 스크립트

| 스크립트 | 용도 | 주요 옵션 |
| --- | --- | --- |
| `convert_myeloid_h5ad_to_seurat.R` | Milo 파이프라인용 **최소** 변환 (counts + 핵심 메타데이터 + 필수 reduction) | `required_meta_cols = c("batch", "Core", "intra_T", "TIL_2", "Annot3")`, `required_obsm_keys = c("X_scVI", "X_umap")` |
| `convert_myeloid_h5ad_with_uns.R` | uns/var 을 유지한 **full** 변환 실험 | `store_var = TRUE`, `store_uns = TRUE`, `required_obsm_keys = NULL` (모든 reduction 유지) |

### 실행 예시
```bash
cd /data/kjc2/git_repo/_wt/h5ad2sobj
Rscript convert_myeloid_h5ad_to_seurat.R          # 최소 버전 (32 MB)
Rscript convert_myeloid_h5ad_with_uns.R          # full 버전 (75 MB)
```

---

## 3. 산출물 비교

| 파일 | 크기 | 구성 | 비고 |
| --- | --- | --- | --- |
| `myeloid_v1_TNBC.qs` | 32 MB | counts, 5개 메타데이터(`batch/Core/intra_T/TIL_2/Annot3`), reductions(`integrated.scvi`, `umap`) | Milo 파이프라인 기본입력 |
| `myeloid_v1_TNBC_full.qs` | 75 MB | 위 + reductions 전체 (`pca`, `spatial` 등), `@misc$anndata_var`(5컬럼), `@misc$anndata_uns`(22 entries) | uns/var 검증용 |

원본 h5ad: 1.1 GB → 최소 변환: 32 MB (약 97% 감소).  
장점: I/O 비용 절감, Milo 파이프라인이 필요한 정보만 적재.

---

## 4. Milo 분석 스크립트

- `run_milo_myeloid_TNBC.R`
  - 입력: `myeloid_v1_TNBC.qs`
  - 비교: `intra_T (low vs high)`, `TIL_2 (low vs high)`
  - 자동으로 low/high 셀만 필터링하고, Seurat reduction을 확인하여 `integrated.scvi` 또는 `pca`를 선택
  - 결과: `milo_results/<comparison>/milo_*` (.qs, plots 등) 및 로그
- `convert_myeloid_h5ad_with_uns.R`로 생성한 full 객체도 동일 스크립트로 사용 가능 (단, 필요시 `seurat_qs_path`만 변경)

---

## 5. 워크트리 및 경로

| 용도 | 경로 |
| --- | --- |
| 메인 repo | `/data/kjc2/git_repo/mylit` |
| h5ad2sobj worktree | `/data/kjc2/git_repo/_wt/h5ad2sobj` |
| milopy worktree | `/data/kjc2/git_repo/_wt/milopy` |

> h5ad2sobj worktree는 기존 `/mylit/_wt/...`에서 `/data/kjc2/git_repo/_wt/...`로 옮겼으며, 스크립트들도 동일 경로를 바라보도록 업데이트했다.

워크트리 이동 예시:
```bash
cd /data/kjc2/git_repo/mylit
git worktree move _wt/h5ad2sobj ../_wt/h5ad2sobj
```

---

## 6. 재현 절차 요약

1. **h5ad 구조 점검 (선택)**
   ```bash
   python inspect_myeloid_h5ad.py
   ```
2. **최소 변환**
   ```bash
   Rscript convert_myeloid_h5ad_to_seurat.R
   ```
3. **full 변환 (필요 시)**
   ```bash
   Rscript convert_myeloid_h5ad_with_uns.R
   ```
4. **Milo 분석**
   ```bash
   Rscript run_milo_myeloid_TNBC.R
   ```

---

## 7. 향후 참고

- `store_obsp`는 h5ad에 obsp가 존재할 때만 동작한다. 현재 TNBC Myeloid 데이터에는 해당 매트릭스가 없으므로 “No obsp keys available to store” 경고만 출력된다.
- milopy worktree(`/_wt/milopy`)는 경로 리팩터링 및 Python 기반 파이프라인 시도(`run_milo_myeloid_TNBC.py`)가 남아있다. 필요 시 메인 저장소에 병합하거나 정리할 것.

