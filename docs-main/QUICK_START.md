# Quick Start: Git Subtree Setup

## 1. GitHub 저장소 생성

1. https://github.com/new 접속
2. Repository name: `mylit-docs`
3. Description: "Documentation repository for mylit project"
4. **Private** 또는 **Public** 선택
5. **중요**: README, .gitignore, license **모두 체크 해제**
6. "Create repository" 클릭

## 2. 저장소 URL 확인

생성 후 SSH URL 복사:
```
git@github.com:JasonK94/mylit-docs.git
```

## 3. 마이그레이션 실행

메인 저장소에서 실행:

```bash
cd /home/user3/data_user3/git_repo/mylit
./scripts/setup-docs-subtree.sh git@github.com:JasonK94/mylit-docs.git
```

스크립트가 자동으로:
- ✅ 현재 `docs-main/` 내용을 새 저장소로 복사
- ✅ 새 저장소에 커밋 및 푸시
- ✅ 메인 저장소에 subtree로 추가
- ✅ Remote 설정

## 4. 검증

```bash
# Subtree remote 확인
git remote -v | grep docs-main

# docs-main/ 파일 확인
git ls-files docs-main/ | head -5

# 동기화 테스트
./scripts/docs-sync.sh pull main
```

## 5. 워크트리 설정 (선택사항)

각 워크트리에서 docs를 사용하려면:

```bash
# 예: fgs 워크트리
cd /data/user3/git_repo/_wt/fgs
git remote add docs-main-origin git@github.com:JasonK94/mylit-docs.git
git fetch docs-main-origin
```

자세한 내용은 `WORKTREE_SETUP.md` 참조.

## 완료!

이제 어느 워크트리에서든 `docs-main/`을 수정하고 동기화할 수 있습니다:

```bash
# 수정 후
git add docs-main/
git commit -m "Update docs"
./scripts/docs-sync.sh push main

# 다른 워크트리에서 최신 버전 가져오기
./scripts/docs-sync.sh pull main
```

