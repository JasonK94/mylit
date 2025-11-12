# 다음 에이전트를 위한 프롬프트

## 현재 상황
R 패키지 `myR`의 gene signature discovery 및 differential expression 분석 함수 개발 중입니다.

## 최근 완료된 작업
1. **`find_gene_signature_v5.2`, `v5.3`**: 완성 및 테스트 완료
2. **`FGS_v5.2`, `FGS_v5.3`**: `signature.R`에 별칭 추가 완료
3. **`runMUSCAT_v5`**: 정상 작동 확인
4. **v5.3 GAM 동적 k**: 정상 작동 확인
5. **pseudobulk 호환성**: v5.2, v5.3 모두 확인

## 현재 이슈

### 1. runMAST_v1 MAST 패키지 호환성 문제 (우선순위 높음)
**문제**: MAST 1.28.0에서 `FromSeurat`가 제거되어 `runMAST_v1`이 작동하지 않음
- 위치: `myR/R/test_analysis.R` (line ~40-105)
- 현재 코드는 `MAST::SceToSingleCellAssay()` 사용하지만 내부적으로 `FromSeurat` 호출
- `FromMatrix` fallback 추가했으나 여전히 오류 발생

**요청사항**:
- MAST 패키지 최신 버전 문서 확인
- `FromMatrix` 직접 사용 또는 다른 변환 방법 구현
- 또는 MAST 대신 다른 DE 방법 제안

### 2. runNEBULA_v1 최종 테스트 (우선순위 중간)
**상태**: 코드 수정 완료 (HL method 추가, fallback 구현)
- 위치: `myR/R/test_analysis.R` (line ~107-275)
- 다운샘플링 데이터로 테스트 필요
- 성공 시 full 데이터로 테스트

## 핵심 파일
- `myR/R/test.R`: `find_gene_signature_v5.2`, `v5.3`
- `myR/R/signature.R`: `FGS_v5.2`, `FGS_v5.3`
- `myR/R/test_analysis.R`: `runMAST_v1`, `runNEBULA_v1`, `runMUSCAT_v5`

## 테스트 데이터
- 다운샘플링: `/data/user3/sobj/IS_scvi_251107_ds2500.qs`
- 원본: `/data/user3/sobj/IS_scvi_251107.qs`

## 작업 환경
- 프로젝트: `/home/user3/data_user3/git_repo/mylit`
- 작업 디렉토리: `/home/user3/GJC_KDW_250721`
- R 세션 초기화: `st/start.R` 자동 실행

## 다음 단계
1. **runMAST_v1 수정**: MAST 패키지 호환성 문제 해결
2. **runNEBULA_v1 테스트**: 다운샘플링 → full 데이터 순서로 테스트
3. **결과 저장**: 모든 테스트 결과는 `qs::qsave()`로 `/data/user3/sobj/`에 저장
4. **커밋**: 작업 완료 후 변경사항 커밋

## 참고
- 모든 함수는 `myR` 패키지의 일부
- `devtools::load_all()`로 패키지 재로드 필요
- 테스트는 다운샘플링 데이터로 먼저 수행 권장

