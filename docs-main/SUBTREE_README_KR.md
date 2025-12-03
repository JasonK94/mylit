# 문서 서브트리 관리

이 디렉토리(`docs-main/`)는 **Git Subtree**로 관리되며, 독립적인 버전 관리를 유지하면서도 모든 워크트리/브랜치에서 편집할 수 있고 병합 충돌을 최소화합니다.

## 구조

```
docs-main/
├── guide.md                    # 메인 프로젝트 가이드 (영어)
├── guide_KR.md                # 메인 프로젝트 가이드 (한국어)
├── README.md                   # 루트 README
├── SUBTREE_README_EN.md       # 이 파일 (영어 버전)
├── SUBTREE_README_KR.md       # 이 파일 (한국어 버전)
├── CONTEXT_DOCUMENTATION.md   # Agent 컨텍스트 문서
├── integrated/
│   ├── INTEGRATED_GUIDE.md     # 종합 모듈 가이드 (영어)
│   └── INTEGRATED_GUIDE_KR.md  # 종합 모듈 가이드 (한국어)
├── rules/
│   └── DOCS_ORGANIZATION_RULE.md  # 문서 조직 규칙
├── dev/
│   ├── DEVLOG.md              # 개발 로그 (영어)
│   ├── DEVLOG_Korean.md       # 개발 로그 (한국어)
│   ├── CHANGELOG.md           # 변경 로그 (영어)
│   ├── CHANGELOG_Korean.md    # 변경 로그 (한국어)
│   ├── context.md             # 프로젝트 컨텍스트 (영어)
│   ├── context_Korean.md      # 프로젝트 컨텍스트 (한국어)
│   └── README_Korean.md       # 패키지 README (한국어)
└── config/
    └── vars_config.R          # 변수 설정 스크립트
```

## 목적

이 디렉토리는 **Git Subtree**로 관리되어:
- **병합 충돌 방지**: 독립적인 저장소로 여러 워크트리가 동시에 편집할 때 충돌을 줄입니다
- **워크트리 편집 가능**: 모든 워크트리/브랜치에서 편집할 수 있습니다
- **독립적인 버전 관리**: 별도 저장소에서 자체 git 히스토리를 가집니다
- **쉬운 동기화**: `git subtree` 명령어를 사용하여 독립적으로 push/pull할 수 있습니다

## 작동 방식

`docs-main/`은 별도 저장소(`mylit-docs`)를 가리키는 Git Subtree입니다:
- **원격 저장소**: `git@github.com:JasonK94/mylit-docs.git` (또는 설정된 원격 저장소)
- **원격 이름**: `docs-main-origin`
- **통합**: 메인 저장소에 포함되어 있지만 독립적으로 동기화할 수 있습니다

## 사용법

### 모든 워크트리에서 편집

```bash
# 모든 워크트리에서
cd /path/to/worktree

# 파일 편집
vim docs-main/guide.md

# 변경사항 커밋 (현재 브랜치에)
git add docs-main/
git commit -m "Update guide.md"

# docs 저장소로 push
./scripts/docs-sync.sh push main
```

### 최신 변경사항 가져오기

```bash
# 원격에서 최신 docs 가져오기
./scripts/docs-sync.sh pull main
```

### Main 브랜치에서만 편집 (더 안전한 방법)

```bash
# main 브랜치로 전환
cd /home/user3/data_user3/git_repo/mylit
git checkout main

# 편집 및 커밋
vim docs-main/guide.md
git add docs-main/
git commit -m "Update guide.md"
git push origin main

# docs 저장소로도 push
./scripts/docs-sync.sh push main
```
### alias 설정

<<<<<<< HEAD
`--prefix` 옵션을 매번 입력하는 대신 git alias를 설정하여 간단하게 사용할 수 있습니다:

```bash
# Alias 설정 (현재 저장소에만)
git config alias.docs-push 'subtree push --prefix=docs-main docs-main-origin main'
git config alias.docs-pull 'subtree pull --prefix=docs-main docs-main-origin main'

# 또는 전역 설정 (모든 저장소에서 사용)
git config --global alias.docs-push 'subtree push --prefix=docs-main docs-main-origin main'
git config --global alias.docs-pull 'subtree pull --prefix=docs-main docs-main-origin main'
```

설정 후 사용:

```bash
# Push
git docs-push

# Pull
git docs-pull
```

=======
>>>>>>> pipe
<<<<<<< HEAD
=======
### alias 설정

`--prefix` 옵션을 매번 입력하는 대신 git alias를 설정하여 간단하게 사용할 수 있습니다:

```bash
# Alias 설정 (현재 저장소에만)
git config alias.docs-push 'subtree push --prefix=docs-main docs-main-origin main'
git config alias.docs-pull 'subtree pull --prefix=docs-main docs-main-origin main'

# 또는 전역 설정 (모든 저장소에서 사용)
git config --global alias.docs-push 'subtree push --prefix=docs-main docs-main-origin main'
git config --global alias.docs-pull 'subtree pull --prefix=docs-main docs-main-origin main'
```

설정 후 사용:

```bash
# Push
git docs-push

# Pull
git docs-pull
```

>>>>>>> 2bdadc73a6b5d3ff4c345f477b3c7b3dcbeb5b9e
## 설정

서브트리를 아직 설정하지 않았다면 다음을 참조하세요:
- `archive/SUBTREE_SETUP_INSTRUCTIONS_EN.md` - 초기 설정 가이드
- `SUBTREE_WORKTREE_SETUP.md` - 워크트리 설정 가이드
- `archive/SUBTREE_STRATEGY_PLAN.md` - 상세 전략 문서

## 충돌 해결

여러 워크트리가 동시에 편집한 경우:

1. 먼저 최신 버전을 pull:
   ```bash
   ./scripts/docs-sync.sh pull main
   ```

2. 필요시 수동으로 충돌 해결

3. 해결된 변경사항 push:
   ```bash
   git add docs-main/
   git commit -m "Resolve conflicts"
   ./scripts/docs-sync.sh push main
   ```

## 참고사항

- push 전에 항상 pull하여 충돌을 피하세요
- 큰 변경보다는 작고 자주 커밋하는 것이 좋습니다
- 서브트리는 전체 git 히스토리를 보존합니다
- 각 워크트리는 `docs-main-origin` 원격이 설정되어 있어야 합니다 (자세한 내용은 `SUBTREE_WORKTREE_SETUP.md` 참조)

