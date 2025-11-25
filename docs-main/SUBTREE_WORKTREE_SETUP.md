# Worktree Setup for Docs Subtree

각 워크트리에서 docs-main subtree를 사용하려면 remote를 추가해야 합니다.

## 자동 설정 (권장)

메인 저장소에서 설정이 완료되면, 각 워크트리에서 다음 명령어를 실행:

```bash
cd /path/to/worktree
git remote add docs-main-origin <docs-repo-url> 2>/dev/null || git remote set-url docs-main-origin <docs-repo-url>
git fetch docs-main-origin
```

## 수동 설정

각 워크트리 디렉토리에서:

```bash
# 예: fgs 워크트리
cd /data/user3/git_repo/_wt/fgs

# Remote 추가 (이미 있으면 skip)
git remote add docs-main-origin git@github.com:JasonK94/mylit-docs.git

# 최신 버전 가져오기
git fetch docs-main-origin

# docs-main 동기화
./scripts/docs-sync.sh pull main
```

## 모든 워크트리 일괄 설정 스크립트

메인 저장소에서 실행:

```bash
cd /home/user3/data_user3/git_repo/mylit

# docs-main-origin URL 가져오기
DOCS_URL=$(git remote get-url docs-main-origin)

# 각 워크트리에 remote 추가
for wt in /data/user3/git_repo/_wt/*; do
    if [ -d "$wt/.git" ]; then
        echo "Setting up: $wt"
        (cd "$wt" && git remote add docs-main-origin "$DOCS_URL" 2>/dev/null || git remote set-url docs-main-origin "$DOCS_URL")
        (cd "$wt" && git fetch docs-main-origin)
    fi
done
```

## 워크트리별 사용법

### 1. 워크트리에서 docs-main 수정

```bash
cd /data/user3/git_repo/_wt/fgs

# docs-main/ 수정
vim docs-main/guide.md

# 변경사항 커밋 (현재 워크트리 브랜치에)
git add docs-main/
git commit -m "Update guide.md from fgs worktree"

# docs-main 저장소로 push
./scripts/docs-sync.sh push main
```

### 2. 다른 워크트리에서 최신 버전 가져오기

```bash
cd /data/user3/git_repo/_wt/lds

# 최신 docs-main 가져오기
./scripts/docs-sync.sh pull main
```

### 3. Main 브랜치에서만 수정 (안전한 방법)

```bash
# Main 브랜치로 전환
cd /home/user3/data_user3/git_repo/mylit
git checkout main

# docs-main/ 수정
vim docs-main/guide.md

# 변경사항 커밋
git add docs-main/
git commit -m "Update guide.md"
git push origin main

# docs-main 저장소로도 push
./scripts/docs-sync.sh push main
```

## 충돌 해결

여러 워크트리에서 동시에 수정한 경우:

```bash
# 1. 최신 버전 pull
./scripts/docs-sync.sh pull main

# 2. 충돌이 있으면 수동으로 해결
# ... 충돌 파일 편집 ...

# 3. 해결 후 다시 push
git add docs-main/
git commit -m "Resolve conflicts in docs-main"
./scripts/docs-sync.sh push main
```

## 주의사항

1. **항상 pull 먼저**: 다른 워크트리에서 수정했을 수 있으므로 push 전에 pull 권장
2. **작은 단위로 commit**: 큰 변경사항은 여러 커밋으로 나누기
3. **충돌 시 신중하게**: 자동 병합이 어려우면 main에서 수동으로 해결

