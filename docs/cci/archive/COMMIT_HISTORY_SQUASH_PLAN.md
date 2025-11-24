# Commit History Squash Plan

## 목적
ligand-target heatmap 버그 수정 과정에서 발생한 여러 commit을 하나의 의미있는 commit으로 축약하여 이력을 정리합니다.

## 현재 상황
버그 수정 과정에서 다음과 같은 commit들이 생성되었습니다:
1. 디버깅 정보 추가
2. 파라미터 이름 수정 (`ligands` → `ligand`)
3. geneset_oi 디버깅 정보 추가
4. replot validation 수정

## 권장 방법

### Option 1: Interactive Rebase (추천)
```bash
# 최근 N개의 commit을 하나로 합치기
git rebase -i HEAD~N

# 또는 특정 commit부터
git rebase -i <base_commit_hash>

# 편집기에서:
# - 첫 번째 commit은 "pick"으로 유지
# - 나머지는 "squash" 또는 "s"로 변경
# - commit 메시지를 통합하여 작성
```

### Option 2: Reset and Recommit
```bash
# 현재 변경사항을 임시 저장
git stash

# 특정 commit으로 reset (soft: 변경사항 유지)
git reset --soft <base_commit_hash>

# 모든 변경사항을 하나의 commit으로
git commit -m "Fix ligand-target heatmap: parameter name bug and replot validation

- Fixed critical bug: changed 'ligands=' to 'ligand=' in get_weighted_ligand_target_links call
  This was causing all ligand-target links to be empty (0 rows)
- Added debugging information for geneset_oi and ligand-target links
- Fixed replot_nichenet_circos validation to accept recordedplot with data in [[2]]
- Improved file saving with dev.flush() before dev.off()
- Heatmap file now generates correctly (62KB vs previous 36KB empty file)"
```

### Option 3: Create New Branch and Cherry-pick
```bash
# 새 브랜치 생성
git checkout -b cci/ligand-target-fix

# 필요한 commit만 cherry-pick
git cherry-pick <commit1> <commit2> <commit3>

# 원래 브랜치로 돌아가서 merge
git checkout cci
git merge --squash cci/ligand-target-fix
git commit -m "통합된 메시지"
```

## 추천 통합 Commit 메시지

```
Fix ligand-target heatmap: parameter name bug and replot validation

Critical bug fix:
- Changed 'ligands=' to 'ligand=' in get_weighted_ligand_target_links call
- This was causing all ligand-target links to be empty (0 rows)
- Heatmap file now generates correctly (62KB vs previous 36KB empty file)

Improvements:
- Added debugging information for geneset_oi and ligand-target links
- Fixed replot_nichenet_circos validation to accept recordedplot with data in [[2]]
- Improved file saving with dev.flush() before dev.off() and longer sleep time

Testing:
- Verified active_ligand_target_links_df now contains data
- Confirmed heatmap file size increased from 36KB to 62KB
```

## 주의사항
1. **이미 push된 commit이면**: `git push --force`가 필요하지만, 협업 중이면 주의 필요
2. **로컬 브랜치만**: `git rebase -i`로 안전하게 수정 가능
3. **백업**: 작업 전 `git branch backup-cci`로 백업 권장

## 실행 시점
- 모든 테스트 완료 후
- replot 기능 확인 후
- 사용자 승인 후

