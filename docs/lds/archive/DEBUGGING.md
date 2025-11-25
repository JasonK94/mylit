# LDS 디버깅 가이드

## 일반적인 오류 및 해결 방법

### 1. 패키지 관련 오류

#### `MulticoreParam` 오류
```
Error: could not find function "MulticoreParam"
```

**해결**:
```r
library(BiocParallel)
# 또는
BiocParallel::MulticoreParam(workers = 4)
```

#### `voomWithDreamWeights` 또는 `dream` 함수를 찾을 수 없음
```
Error: could not find function "voomWithDreamWeights"
```

**해결**:
```r
# limma 패키지가 최신 버전인지 확인
packageVersion("limma")
# 최신 버전 설치
BiocManager::install("limma")
library(limma)
```

#### `sva` 함수를 찾을 수 없음
```
Error: could not find function "sva"
```

**해결**:
```r
BiocManager::install("sva")
library(sva)
```

### 2. 데이터 형식 오류

#### Seurat 객체에서 count 추출 실패
```
Error in GetAssayData(sobj, layer = layer)
```

**해결**:
```r
# layer 확인
DefaultAssay(sobj)
# count 데이터 확인
GetAssayData(sobj, layer = "counts")[1:5, 1:5]
# 다른 layer 시도
GetAssayData(sobj, layer = "RNA")[1:5, 1:5]
```

#### 메타데이터 불일치
```
Error: Count matrix의 열 수(샘플)와 meta.data의 행 수(샘플)가 일치하지 않습니다.
```

**해결**:
```r
# 확인
ncol(GetAssayData(sobj, layer = "counts"))
nrow(sobj@meta.data)

# 일치 확인
all(colnames(GetAssayData(sobj, layer = "counts")) == rownames(sobj@meta.data))
```

### 3. 포뮬러 관련 오류

#### lme4 포뮬러 파싱 실패
```
Error in lme4::nobars(formula)
```

**해결**:
```r
# 포뮬러 확인
formula <- ~ g3_clean + (1|hos_no) + (1|GEM)
class(formula)  # "formula"여야 함

# 변수 확인
all(c("g3_clean", "hos_no", "GEM") %in% colnames(sobj@meta.data))
```

#### 변수가 메타데이터에 없음
```
Error: variable 'hos_no' not found
```

**해결**:
```r
# 메타데이터 확인
colnames(sobj@meta.data)
# 변수명 확인 (대소문자 주의)
```

### 4. SVA 관련 오류

#### SVA가 SV를 찾지 못함
```
... SVA가 유의미한 대리 변수(SV)를 찾지 못했습니다.
```

**원인**:
- 샘플 수가 너무 적음
- 고정 효과 모델이 이미 대부분의 변동성을 설명함
- 정상적인 경우일 수 있음

**해결**:
- SV 없이 진행됨 (정상)
- 샘플 수를 늘리거나 다른 데이터 사용

#### SVA 실행 중 오류
```
Error in sva(v_sva$E, mod = mod_sva, mod0 = mod0_sva, n.sv = NULL)
```

**해결**:
```r
# mod와 mod0 확인
mod_sva <- model.matrix(~ g3_clean, data = meta.data)
mod0_sva <- model.matrix(~ 1, data = meta.data)

# 차원 확인
dim(mod_sva)
dim(mod0_sva)
dim(v_sva$E)

# 샘플 수 확인 (너무 적으면 문제)
nrow(meta.data)
```

### 5. dream 관련 오류

#### dream 피팅 실패
```
Error in dream(v_dream, final_formula, meta.data, BPPARAM = BPPARAM_SETUP)
```

**해결**:
```r
# 포뮬러 확인
final_formula
# 메타데이터에 모든 변수가 있는지 확인
all.vars(final_formula) %in% colnames(meta.data)

# 샘플 수 확인 (너무 적으면 문제)
nrow(meta.data)
```

#### 메모리 부족
```
Error: cannot allocate vector of size X Mb
```

**해결**:
- 더 작은 데이터셋 사용
- 유전자 필터링 강화
- `n_cores` 줄이기

### 6. 디버깅 팁

#### 단계별 실행
```r
# 1. 데이터 확인
head(is5s_test@meta.data[, c("g3_clean", "hos_no", "GEM")])

# 2. 포뮬러 테스트
formula <- ~ g3_clean + (1|hos_no) + (1|GEM)
lme4::nobars(formula)

# 3. DGEList 생성 테스트
counts <- GetAssayData(is5s_test, layer = "counts")
dge <- DGEList(counts, samples = is5s_test@meta.data)

# 4. SVA 테스트
mod <- model.matrix(~ g3_clean, data = is5s_test@meta.data)
mod0 <- model.matrix(~ 1, data = is5s_test@meta.data)
v <- voom(dge, mod, plot = FALSE)
sva_result <- sva(v$E, mod = mod, mod0 = mod0, n.sv = NULL)
```

#### 로그 확인
LDS 함수는 각 단계마다 메시지를 출력합니다:
```
1/7: 입력 데이터 처리 중...
2/7: 포뮬러 파싱 중...
3/7: DGEList 생성 및 필터링 중...
4/7: SVA 실행 (숨겨진 변동성 탐색)...
5/7: 최종 모델 포뮬러 생성 중...
6/7: limma-dream 실행 (Core: X개)...
7/7: 분석 완료.
```

어느 단계에서 오류가 발생하는지 확인하세요.

## 문제 리포트

오류가 발생하면 다음 정보를 수집하세요:

1. **오류 메시지 전체**
2. **데이터 정보**:
   - 셀 수, 유전자 수
   - 샘플 수 (hos_no)
   - 배치 수 (GEM)
   - g3_clean 분포
3. **실행 환경**:
   - R 버전
   - 패키지 버전 (limma, edgeR, lme4, sva)
4. **어느 단계에서 실패했는지** (LDS 함수의 1/7 ~ 7/7 중)

