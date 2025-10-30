# 작업 요약 (Work Summary)

## 1. 프로젝트 목표

`myR/function_analysis.md` 마크다운 파일에 정리된 R 함수 목록을 분석하여, 컴퓨터가 처리하기 용이한 CSV 형태의 데이터로 변환하고, 이를 한국어로 번역하는 데이터 파이프라인을 구축합니다.

## 2. 작업 절차 및 산출물

### 2.1. 마크다운 파싱 (`01_parse_md_to_csv.py`)

- **기능**: `myR/function_analysis.md` 파일을 읽어들여 구조를 분석하고, 각 함수에 대한 정보를 추출합니다.
- **실행**: `python mypy/package_tracking/01_parse_md_to_csv.py`
- **산출물**: `function_list_en.csv`
  - `id`, `directory`, `filepath`, `function_name`, `description`, `input_parameters`, `returns` 컬럼을 포함하는 영문 데이터 파일입니다.

### 2.2. CSV 번역 (`02_translate_csv_with_deepl.py`)

- **기능**: 생성된 `function_list_en.csv` 파일의 주요 텍스트 컬럼들(`description`, `input_parameters`, `returns`)을 DeepL API를 사용하여 한국어로 번역합니다.
- **실행**: `python mypy/package_tracking/02_translate_csv_with_deepl.py`
- **산출물**: `function_list_ko.csv`
  - 영문 CSV와 동일한 구조를 가지며, 주요 내용이 한국어로 번역된 파일입니다.
  - **참고**: DeepL 무료 API의 월간 사용량 제한으로 인해, 현재 파일은 일부만 번역된 상태일 수 있습니다.

## 3. 스크립트 및 환경 설정

- **스크립트 위치**: 모든 파이썬 스크립트는 `mypy/package_tracking/` 디렉터리 내에 위치합니다.
- **API 키 관리**:
  - DeepL API 키는 프로젝트 루트의 `.env` 파일에 `DEEPL_API_KEY=...` 형식으로 저장됩니다.
  - 보안을 위해 `.env` 파일은 `.gitignore`에 등록되어 Git 버전 관리에서 제외됩니다.
- **필요 라이브러리**:
  - `deepl`: DeepL API 사용
  - `python-dotenv`: `.env` 파일 로드
  - 설치: `pip install deepl python-dotenv`
