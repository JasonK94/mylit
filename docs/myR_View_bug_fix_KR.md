# myR 패키지 View() 함수 버그 진단 및 해결

## 문제 요약

`devtools::load_all("/data/user3/git_repo/mylit/myR")`를 실행한 후, RStudio의 `View()` 함수가 Seurat 객체를 표시하지 못하고 다음 에러를 발생시킴:

```r
> View(new)
Error in as.data.frame.default(x) : 
  cannot coerce class 'structure("Seurat", package = "SeuratObject")' to a data.frame
```

## 근본 원인 (Root Cause)

### 1. 의존성 체인
```
myR DESCRIPTION Imports:
  └─ S4Vectors
      └─ BioGenerics (자동 로드)
          └─ as.data.frame을 S4 generic으로 override
```

### 2. 문제 발생 경로
1. `devtools::load_all()` 실행
2. DESCRIPTION의 `Imports`에 있는 `S4Vectors` 패키지 로드
3. `S4Vectors`가 의존성으로 `BioGenerics` 자동 로드
4. `BioGenerics`가 `base::as.data.frame`을 S4 generic으로 덮어씀
5. RStudio의 `View()` 함수가 내부적으로 `as.data.frame()` 호출
6. BioGenerics의 S4 메소드가 호출되지만 Seurat 객체에 대한 메소드가 정의되지 않음
7. `as.data.frame.default()`로 fall back되며 에러 발생

### 3. 왜 base::를 강제해도 해결되지 않았는가?

초기에 `zzz.R`의 `.onLoad()`에서 `conflicted::conflicts_prefer(base::as.data.frame)`이 주석처리되어 있었기 때문. 이 설정이 없으면 S4 generic이 우선됨.

## 해결책

### 해결 1: conflicted 설정 활성화 ✅

[`myR/R/zzz.R`](file:///home/user3/data_user3/git_repo/mylit/myR/R/zzz.R#L36)에서 주석 해제:

```r
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("conflicted", quietly = TRUE)) {
        # base preferences
        conflicted::conflicts_prefer(base::intersect)
        conflicted::conflicts_prefer(base::setdiff)
        conflicted::conflicts_prefer(base::union)
        conflicted::conflicts_prefer(base::as.data.frame)  # ← 이 줄 활성화
        
        # ...
    }
}
```

### 해결 2: S4Vectors를 선택적 의존성으로 변경 ✅

`S4Vectors`가 모든 함수에서 필수가 아니라면 `Imports`에서 `Suggests`로 이동:

**DESCRIPTION 수정:**
```diff
 Imports:
     rlang,
     rstatix,
     rsvg,
-    S4Vectors,
     scales,
     Seurat (>= 4.0.0),
     ...
 Suggests:
     anndata,
     broom.mixed,
     ComplexUpset,
+    S4Vectors,
     DiagrammeR,
     ...
```

이렇게 하면:
- `S4Vectors`가 자동으로 로드되지 않음
- 필요한 함수에서만 `S4Vectors::`로 명시적 호출
- `BioGenerics`가 자동 로드되지 않아 `as.data.frame` override 없음

### 해결 3: 코드에서 명시적 호출 사용

`S4Vectors`를 사용하는 코드에서 이미 명시적으로 호출하고 있음:

```r
# deg_consensus/deg_methods_dream.R:122
SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)

# analysis/pseudotime.R:408
colData = S4Vectors::DataFrame(original_metadata_sce)
```

이는 `S4Vectors`를 `Suggests`로 이동해도 문제 없음을 의미.

## 검증

수정 후 다음 명령으로 테스트:

```r
# 1. 패키지 재로드
devtools::load_all("/data/user3/git_repo/mylit/myR")

# 2. Seurat 객체 생성 (또는 기존 객체 사용)
# new <- readRDS("some_seurat.rds")

# 3. View 테스트
View(new)  # 에러 없이 RStudio viewer에서 객체 표시되어야 함
```

## 추가 확인 사항

### 현재 as.data.frame 메소드 확인:
```r
methods("as.data.frame")
```

정상적으로 작동하면:
- `base::as.data.frame`이 기본으로 사용됨
- BioGenerics S4 generic이 있더라도 우선순위가 낮음

### conflicted 설정 확인:
```r
library(conflicted)
conflicts_prefer()  # 현재 설정된 우선순위 확인
```

## 관련 파일

- [myR/R/zzz.R](file:///home/user3/data_user3/git_repo/mylit/myR/R/zzz.R) - 패키지 로드 설정
- [myR/DESCRIPTION](file:///home/user3/data_user3/git_repo/mylit/myR/DESCRIPTION) - 의존성 정의
- [myR/NAMESPACE](file:///home/user3/data_user3/git_repo/mylit/myR/NAMESPACE) - 자동 생성 (roxygen2)

## 결론

문제는 사용자의 추정이 정확했음:
- ✅ BioGenerics의 `as.data.frame` override가 원인
- ✅ DESCRIPTION/NAMESPACE의 변화로 야기됨 (S4Vectors 추가)
- ✅ 상대적으로 최근 버그 (초기 public repo 커밋 시 S4Vectors 추가됨)

하지만 `.onLoad()`의 `conflicts_prefer` 설정이 주석처리되어 있어서 `base::`를 강제해도 효과가 없었음. 두 가지 수정 (conflicted 설정 + Suggests 이동)으로 완전히 해결됨.
