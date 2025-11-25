# Documentation Subtree Strategy Plan

## 요구사항 정리

1. **문제**: 여러 워크트리에서 `docs-main/`을 동시 수정 시 merge conflict 빈발
2. **목표**: 어느 브랜치/워크트리에서 시작했든 main에 적용하고 다른 워크트리에도 동기화
3. **제약**: 별도 폴더만으로는 충돌 해결 안 됨 → 별도 Git 저장소 필요
4. **선호**: 각 워크트리에서 독립 수정 + 즉시 commit & push로 동기화

## 전략 비교

### Option 1: Git Subtree (추천 ⭐)

**구조:**
- `docs-main/`을 별도 Git 저장소로 분리
- 메인 저장소에 `git subtree`로 통합
- 각 워크트리에서 `git subtree pull/push`로 동기화

**장점:**
- ✅ 충돌 최소화: 별도 저장소에서 독립 버전 관리
- ✅ 워크트리별 독립 수정 가능
- ✅ 즉시 push/pull로 동기화 가능
- ✅ Git 히스토리 보존
- ✅ 메인 저장소에 통합되어 있어 별도 clone 불필요

**단점:**
- ⚠️ `git subtree` 명령어 학습 필요
- ⚠️ 충돌 시 수동 해결 필요 (하지만 빈도는 크게 감소)

**워크플로우:**
```bash
# 워크트리에서 수정 후
cd /path/to/worktree
git subtree push --prefix=docs-main docs-main-origin main

# 다른 워크트리에서 최신 버전 가져오기
cd /path/to/other-worktree
git subtree pull --prefix=docs-main docs-main-origin main
```

### Option 2: External Repo + Symlink

**구조:**
- `docs-main/`을 완전히 별도 저장소로 분리
- 메인 저장소 밖에 위치 (예: `/data/user3/git_repo/mylit-docs`)
- 메인 저장소에서 심볼릭 링크로 참조

**장점:**
- ✅ 완전히 독립적인 버전 관리
- ✅ 충돌 없음 (메인 저장소와 완전 분리)
- ✅ 표준 Git 워크플로우 사용 가능

**단점:**
- ❌ 심볼릭 링크 관리 복잡
- ❌ 각 워크트리마다 별도 clone 필요
- ❌ 메인 저장소와 분리되어 있어 통합성 낮음

**워크플로우:**
```bash
# 각 워크트리에서 별도 docs 저장소 clone
cd /path/to/worktree
cd ../
git clone <docs-repo-url> mylit-docs
cd worktree
ln -s ../mylit-docs docs-main

# 수정 후
cd ../mylit-docs
git add . && git commit -m "..." && git push
```

### Option 3: Main-Only Editing

**구조:**
- `docs-main/`을 메인 저장소 내 일반 디렉토리로 유지
- main 브랜치에서만 수정 허용
- 워크트리에서는 읽기 전용

**장점:**
- ✅ 가장 단순한 구조
- ✅ 충돌 없음 (main에서만 수정)

**단점:**
- ❌ 워크트리에서 직접 수정 불가
- ❌ 문제 인식 → main 전환 → 수정 → 워크트리로 복귀 (비효율적)

## 추천: Git Subtree (Option 1)

**이유:**
1. 충돌 최소화하면서도 워크트리에서 직접 수정 가능
2. 메인 저장소에 통합되어 있어 관리 용이
3. 즉시 동기화 가능
4. Git 히스토리 보존

## 구현 계획

### Phase 1: 별도 저장소 생성 및 초기화

1. **새 저장소 생성**
   - GitHub에 `mylit-docs` 저장소 생성 (또는 기존 저장소 사용)
   - 또는 로컬에 bare 저장소 생성

2. **현재 docs-main/ 내용을 새 저장소로 마이그레이션**
   ```bash
   # 임시로 docs-main 내용을 새 저장소로 복사
   cd /tmp
   git clone <docs-repo-url> mylit-docs-temp
   cp -r /home/user3/data_user3/git_repo/mylit/docs-main/* mylit-docs-temp/
   cd mylit-docs-temp
   git add .
   git commit -m "Initial docs-main migration"
   git push
   ```

### Phase 2: 메인 저장소에 Subtree 추가

1. **메인 저장소에서 기존 docs-main/ 제거**
   ```bash
   cd /home/user3/data_user3/git_repo/mylit
   git rm -r docs-main/
   git commit -m "Remove docs-main/ before subtree migration"
   ```

2. **Subtree로 추가**
   ```bash
   git subtree add --prefix=docs-main docs-main-origin main --squash
   # 또는 히스토리 보존
   git subtree add --prefix=docs-main docs-main-origin main
   ```

3. **Remote 추가**
   ```bash
   git remote add docs-main-origin <docs-repo-url>
   ```

### Phase 3: 워크트리 설정

각 워크트리에서:
```bash
cd /path/to/worktree
git remote add docs-main-origin <docs-repo-url>  # 이미 있으면 skip
git fetch docs-main-origin
```

### Phase 4: 워크플로우 스크립트 생성

**Helper 스크립트 생성: `scripts/docs-sync.sh`**
```bash
#!/bin/bash
# docs-main 동기화 스크립트

ACTION=${1:-pull}  # pull or push
BRANCH=${2:-main}

if [ "$ACTION" = "push" ]; then
    echo "Pushing docs-main changes to remote..."
    git subtree push --prefix=docs-main docs-main-origin $BRANCH
elif [ "$ACTION" = "pull" ]; then
    echo "Pulling latest docs-main from remote..."
    git subtree pull --prefix=docs-main docs-main-origin $BRANCH --squash
else
    echo "Usage: $0 [push|pull] [branch]"
    exit 1
fi
```

### Phase 5: 문서화 업데이트

1. `docs-main/README_SUBTREE.md` 업데이트
2. `CONTEXT_DOCUMENTATION.md` 생성 (새 전략 설명)
3. 워크플로우 가이드 작성

## 사용 시나리오

### 시나리오 1: 워크트리에서 직접 수정

```bash
# fgs 워크트리에서 작업 중
cd /data/user3/git_repo/_wt/fgs

# docs-main/ 수정
vim docs-main/guide.md

# 변경사항 커밋 (워크트리 브랜치에)
git add docs-main/
git commit -m "Update guide.md from fgs worktree"

# docs-main 저장소로 push
./scripts/docs-sync.sh push main

# 다른 워크트리에서 최신 버전 가져오기
cd /data/user3/git_repo/_wt/lds
./scripts/docs-sync.sh pull main
```

### 시나리오 2: Main에서만 수정 (안전한 방법)

```bash
# main 브랜치로 전환
cd /home/user3/data_user3/git_repo/mylit
git checkout main

# docs-main/ 수정
vim docs-main/guide.md

# 변경사항 커밋 및 push
git add docs-main/
git commit -m "Update guide.md"
git push origin main

# docs-main 저장소로도 push
./scripts/docs-sync.sh push main

# 워크트리에서 최신 버전 가져오기
cd /data/user3/git_repo/_wt/fgs
./scripts/docs-sync.sh pull main
```

## 충돌 해결 전략

여러 워크트리에서 동시에 수정할 경우:

1. **먼저 push한 변경사항이 우선**
2. **나중에 push하려는 워크트리에서:**
   ```bash
   # 최신 버전 pull
   ./scripts/docs-sync.sh pull main
   
   # 충돌 해결
   # ... 수동으로 충돌 해결 ...
   
   # 다시 push
   ./scripts/docs-sync.sh push main
   ```

## 다음 단계

1. 사용자 확인: 별도 저장소 URL 또는 생성 방법
2. 구현 시작: Phase 1부터 순차 진행
3. 테스트: 한 워크트리에서 수정 → push → 다른 워크트리에서 pull 테스트

