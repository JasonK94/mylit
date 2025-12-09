# Documentation Organization Summary

## Completed Actions (2025-12-09)

파이프라인 문서 정리를 완료하였습니다. 다음과 같이 구성되었습니다:

### 1. 이동된 파일들

#### Root → docs/pipe/
- `context.md` → `docs/pipe/CONTEXT.md` (파이프라인 현황 및 이슈)

#### Root → docs/pipe/archive/
- `install_log.txt`, `install_log_2.txt`, `install_log_3.txt` (설치 로그)
- `run_pipeline_*.log` (실행 로그 - full, hto, test 등)
- `run_pipeline_*.sh` (실행 스크립트)
- `run_test_v18.log` (테스트 로그)

#### scripts/ → docs/pipe/archive/old_scripts/
- `prep_metadata.R` - 메타데이터 전처리 (메인 파이프라인에 통합됨)
- `fix_manifest_hto.R` - HTO manifest 수정 스크립트
- `update_manifest.R` - Manifest 업데이트 유틸리티
- `update_config_columns.R` - 설정 컬럼 동기화
- `verify_meta.R` - 메타데이터 검증
- `create_downsampled_data.R` - 테스트용 다운샘플링
- `vars_config.R` - 설정 변수 관리
- `load_renv.R` - renv 로딩
- `start_pipe.R` - 구버전 시작 스크립트 (scripts/pipe/start_pipe.R로 대체됨)

### 2. 생성된 문서들

#### docs/pipe/README.md
- 모든 파이프라인 스크립트와 문서의 인덱스
- Core pipeline scripts 위치 안내
- Utility scripts 위치 안내
- 문서 간 상호 참조 링크 정리
- Configuration examples 안내

#### docs/pipe/archive/README.md
- Archive 폴더 내용 설명
- Old scripts의 용도와 현재 대체 방법 안내
- 보관 이유 명시

### 3. 업데이트된 문서

#### docs/pipe/CONTEXT.md
- 상대 경로 링크로 변경 (예: `docs/pipe/COMMANDS.md` → `[COMMANDS.md](COMMANDS.md)`)
- README.md 링크 추가

### 4. 현재 파일 구조

```
pipe/
├── .Rprofile
├── .gitignore
├── config/                    # 설정 파일들
│   ├── manifest_*.csv
│   ├── run_config_*.json
│   └── meta_data_*.csv
├── docs/
│   ├── pipe/                  # 파이프라인 문서 (메인)
│   │   ├── README.md         # **새로 생성** - 전체 인덱스
│   │   ├── CONTEXT.md        # **이동** - 현황 및 이슈
│   │   ├── PIPE_INTEGRATED_GUIDE_KR.md  # 통합 가이드
│   │   ├── COMMANDS.md
│   │   ├── DATA_FLOW.md
│   │   └── archive/          # 보관 문서
│   │       ├── README.md     # **새로 생성**
│   │       ├── old_scripts/  # **새로 생성** - 구버전 스크립트들
│   │       ├── *.log         # **이동** - 실행 로그들
│   │       └── *.txt         # **이동** - 설치 로그들
│   ├── fgs/
│   ├── cci/
│   └── ...
├── scripts/
│   ├── pipe/                 # 현재 파이프라인 실행 스크립트
│   │   ├── pipe1_read_demulti.R
│   │   ├── pipe2_nmz_clustering.R
│   │   └── ...
│   └── pipe_wrapper.sh
├── myR/                      # R 유틸리티 패키지
│   └── R/
│       ├── pipe_utils.R
│       └── pipe_demulti.R
└── logs/                     # 실행 로그 출력 디렉토리
```

### 5. 상호 참조 확인

모든 문서 내 상대 링크를 확인하였으며, 다음과 같이 정리되었습니다:
- `docs/pipe/` 내의 문서들은 상대 경로 사용
- README.md에서 모든 주요 문서로 링크 제공
- Cross-reference가 깨지지 않도록 확인 완료

### 6. 보존된 중요 파일

다음 파일들은 루트에 유지됩니다 (활발히 사용 중):
- `.Rprofile` - R 환경 설정
- `.gitignore` - Git 설정
- `config/` - 현재 사용 중인 설정 파일
- `scripts/pipe/` - 현재 파이프라인 실행 스크립트
- `myR/` - 활발히 개발 중인 R 패키지

## 메인 브랜치 병합 시 주의사항

1. **Archive 폴더**: 오래된 로그와 스크립트는 archive에 보관되어 있으므로, 필요시 참조 가능
2. **문서 링크**: 모든 문서 간 링크가 상대 경로로 수정되어 있어 구조 변경에도 안정적
3. **README.md**: 파이프라인 문서의 진입점으로 활용 가능

---
**정리 완료일**: 2025-12-09
