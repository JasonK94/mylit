# Documentation Organization Rules

이 문서는 프로젝트 내의 문서들을 체계적으로 관리하고 정리하기 위한 규칙을 정의합니다.

## 1. 기본 원칙 (Principles)
1.  **Single Source of Truth**: 모듈별로 핵심 내용을 담은 하나의 **통합 가이드(Integrated Guide)**를 유지한다.
2.  **Separation of Concerns**: 사용자 가이드(How-to)와 개발/디버깅 로그(Log)를 분리한다.
3.  **Archiving**: 해결된 이슈나 오래된 계획 문서는 삭제하지 않고 `archive` 폴더로 이동하여 보존한다.
4.  **Localization**: 가능하면 국문 가이드를 병행하여 작성한다.

## 2. 폴더 구조 (Directory Structure)

```text
docs/
├── [module_name]/                  # 예: fgs, pt.umap
│   ├── [MODULE]_INTEGRATED_GUIDE.md    # 영문 통합 가이드 (Main)
│   ├── [MODULE]_INTEGRATED_GUIDE_KR.md # 국문 통합 가이드
│   ├── README.md                       # 가이드 링크 및 인덱스
│   └── archive/                        # 과거 문서 보관소
│       ├── BUGFIX_...md
│       ├── PLAN_...md
│       └── ...
└── DOCS_ORGANIZATION_RULE.md       # 이 파일
```

## 3. 통합 가이드 작성 가이드 (Integrated Guide Template)

통합 가이드는 다음 섹션을 포함해야 합니다.

### 1. Introduction (소개)
*   모듈의 목적과 핵심 기능을 요약.
*   주요 용어 정의.

### 2. Workflow Visualization (시각화)
*   Mermaid를 사용한 파이프라인/로직 시각화.
*   **중요: 특수 문자 처리 규칙**:
    *   **파이프 문자 (`|`)**: `1|patient` → `["1|patient"]` 또는 `["1 pipe patient"]`로 변경
    *   **등호 (`=`)**: 수식 표현 시 따옴표로 감싸기 `["A_g = mean(S_gm)"]`
    *   **언더스코어 (`_`)**: 변수명에 포함된 경우 따옴표 권장 `["A_g"]`
    *   **기타 특수문자**: `()`, `,`, `/` 등이 포함된 텍스트는 반드시 `["Text"]` 형식으로 작성
    *   **예시**: `["Random Effect<br/>1|GEM/hos_no"]` → `["Random Effect<br/>1 pipe GEM slash hos_no"]` 또는 `["Random Effect<br/>1|GEM/hos_no"]` (따옴표로 감싸기)
*   스타일 가이드:
    *   `flowchart TD` (Top-Down) 권장.
    *   `subgraph`를 사용하여 논리적 블록 구분.
    *   노드 텍스트에 특수 문자가 포함되면 반드시 `["..."]` 형식 사용.

### 3. Development Log & Improvements (개발 로그)
*   주요 버전별 변경 사항 요약.
*   해결된 주요 버그와 해당 커밋/PR 정보.

### 4. User Guide & Warnings (사용자 가이드)
*   **Critical Warnings**: CPU 제한(`taskset`), 메모리 요구사항, 환경 설정(`renv`) 등 필수 주의사항.
*   **Usage**: 핵심 함수 사용 예제 코드.

### 5. Methodology (방법론)
*   알고리즘 및 로직에 대한 상세 설명.

### 6. Appendix (부록)
*   메서드 목록, 파라미터 테이블 등 참조 정보.

## 4. 아카이빙 규칙 (Archiving Rules)

다음과 같은 문서는 `archive/` 폴더로 이동합니다.
*   특정 버그 수정을 위한 일회성 분석 로그 (`BUGFIX_*.md`)
*   이미 구현이 완료된 개발 계획서 (`PLAN_*.md`)
*   진행 상황 보고서 (`PROGRESS_*.md`)
*   특정 시점의 성능 테스트 결과

## 5. README.md 관리
*   해당 모듈의 `README.md`는 통합 가이드로 연결되는 링크를 최상단에 제공해야 합니다.
*   아카이브 폴더에 대한 안내를 포함합니다.

## 6. 문서 Subtree 관리 (Documentation Subtree Management)

### 6.1 핵심 문서 파일 관리
다음 파일들은 `main` 브랜치에서만 관리되며, `docs-main/` 디렉터리에 저장됩니다:

**루트 레벨:**
- `guide.md`, `guide_KR.md`, `README.md`

**docs/ 레벨:**
- `docs/INTEGRATED_GUIDE.md`, `docs/INTEGRATED_GUIDE_KR.md`, `docs/DOCS_ORGANIZATION_RULE.md`

**myR/ 레벨:**
- `myR/DEVLOG.md`, `myR/DEVLOG_Korean.md`, `myR/CHANGELOG.md`, `myR/CHANGELOG_Korean.md`
- `myR/context.md`, `myR/context_Korean.md`, `myR/README_Korean.md`

**scripts/ 레벨:**
- `scripts/vars_config.R`

### 6.2 Subtree 구조
```
docs-main/
├── guide.md
├── guide_KR.md
├── README.md
├── integrated/
│   ├── INTEGRATED_GUIDE.md
│   └── INTEGRATED_GUIDE_KR.md
├── rules/
│   └── DOCS_ORGANIZATION_RULE.md
├── dev/
│   ├── DEVLOG.md
│   ├── DEVLOG_Korean.md
│   ├── CHANGELOG.md
│   ├── CHANGELOG_Korean.md
│   ├── context.md
│   ├── context_Korean.md
│   └── README_Korean.md
└── config/
    └── vars_config.R
```

### 6.3 규칙
1. **main 브랜치 전용**: 이 파일들은 `main` 브랜치에서만 수정 가능합니다.
2. **Symlink 사용**: 원래 위치의 파일들은 `docs-main/`로의 symlink입니다.
3. **보호 메커니즘**: 
   - `.gitattributes`로 파일 속성 지정
   - `pre-commit` hook으로 비-main 브랜치에서의 수정 방지
4. **Feature 브랜치**: Feature 브랜치에서는 이 파일들을 수정하지 않습니다.

### 6.4 Symlink 생성
Symlink는 `create_symlinks.sh` 스크립트를 실행하여 생성합니다:
```bash
bash create_symlinks.sh
```

### 6.5 파일 수정 절차
1. `main` 브랜치로 전환: `git checkout main`
2. `docs-main/` 내 파일 직접 수정 또는 symlink를 통한 수정
3. 변경사항 커밋 및 푸시
4. Feature 브랜치에서는 이 파일들을 수정하지 않음

