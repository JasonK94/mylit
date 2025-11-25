# 변경 이력

`myR`의 모든 주요 변경 사항을 기록합니다.  
- 형식: [Keep a Changelog](https://keepachangelog.com/ko/1.0.0/)  
- 버전 정책: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [발표 예정]

### 추가
- 시그니처 워크플로우 v5.2와 검증 유틸리티, 문서 갱신 (`2d9c50e`, `39a6159`).
- Milo 차등 풍부도 파이프라인(`milo_opus6`)과 메타 수준 유전자 중요도 분석 (`18f1ad1`, `b260bb1`, `a360837`, `7bf2904`).
- TestLISI, PERMANOVA, PCP, PTMFL 등 분석 테스트 및 메타 러너 병렬 기본값 개선 (`25b205b`, `5863176`).
- 의사벌크 DEG 모듈 `pb_deg.R`와 마이그레이션 가이드 (Claude 브랜치 반영, `881589c` 계열).
- 프로젝트 부트스트랩 문서(`context.md`, `NEXT_STEPS.md`)와 다국어 기록 추가 (`cf6d13e`, `18c2810`).

### 변경
- Seurat `slot` → `layer` 전환으로 v5와 호환성 확보 (`881589c`, `39a6159`).
- R 스크립트를 재구성하고 폐기 예정 모듈을 격리, MUSCAT/NEBULA 유틸 표준화 (`497eab8`, `bc19660`).
- `.gitignore`를 강화해 대형 산출물이 저장소에 포함되지 않도록 조정 (`80f7941`, `9774761`).

### 제거
- `projects/`, `trash/` 등 대용량 산출물을 버전 관리에서 제외 (`80f7941`).
- 검증되지 않은 MUSCAT/MAST 실험 코드를 정리하고 안정된 루틴만 유지 (`80260af`).

### 수정
- DESCRIPTION에 `NMF` 의존성을 명시해 `nmf()` 사용 시 오류를 해결 (`e785399`).
- PERMANOVA 메타데이터 처리와 rownames 안정성을 개선 (`9cb9d3a`).


