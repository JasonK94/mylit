# ============================================================================
# CCI 분석 실행 시간 추정
# ============================================================================

# 테스트 결과 기반 추정
# - 20개 DEG, 20개 ligands, 200 targets/ligand (기본값)
#   → Ligand activity: 3.2초
#   → Ligand-target inference: 거의 0초 (0 links found였지만 실제로는 계산됨)

# 사용자 시나리오: 4637개 DEG, 20개 ligands, 50 targets/ligand

n_degs_small <- 20
n_degs_large <- 4637
n_ligands <- 20
top_targets_small <- 200
top_targets_large <- 50

# 1. Ligand Activity Prediction 시간 추정
# predict_ligand_activities는 geneset_oi (DEG) 크기에 대략 비례
# 하지만 실제로는 ligand_target_matrix 크기가 더 중요할 수 있음
ligand_activity_time_small <- 3.2  # 초
ligand_activity_time_large <- ligand_activity_time_small * (n_degs_large / n_degs_small)
cat("=== 시간 추정 ===\n")
cat("1. Ligand Activity Prediction:\n")
cat("   - 20 DEGs: ", ligand_activity_time_small, "초\n")
cat("   - 4637 DEGs: ", round(ligand_activity_time_large, 1), "초 (", round(ligand_activity_time_large/60, 1), "분)\n")
cat("   - 추정 근거: DEG 수에 비례 (선형 스케일링 가정)\n\n")

# 2. Ligand-Target Inference 시간 추정
# get_weighted_ligand_target_links는:
# - 각 ligand마다 실행 (20개)
# - geneset_oi 크기에 비례하지만, 실제로는 ligand_target_matrix에서 필터링하므로 sub-linear
# - top_n_targets_per_ligand에 비례

# 코드의 휴리스틱: ligands × targets / 800 = 분
heuristic_minutes <- (n_ligands * top_targets_large) / 800
cat("2. Ligand-Target Inference (코드 휴리스틱):\n")
cat("   - 20 ligands × 50 targets / 800 = ", round(heuristic_minutes, 1), "분\n")
cat("   - 실제 시간: ", round(heuristic_minutes * 60, 0), "초\n\n")

# 하지만 DEG가 많으면 더 오래 걸릴 수 있음
# 실제 테스트에서 0초였던 것은 links가 0개여서였을 가능성
# 실제로는 각 ligand마다 geneset_oi를 ligand_target_matrix와 매칭하는 작업이 필요

# 보수적 추정: DEG 수에 따라 추가 시간
# geneset_oi가 클수록 매칭 시간이 증가하지만, ligand_target_matrix는 고정 크기
# 실제로는 sub-linear scaling (로그 스케일)
deg_ratio <- n_degs_large / n_degs_small
# 보수적으로 sqrt 스케일링 가정 (sub-linear)
scaling_factor <- sqrt(deg_ratio)  # 약 15.2배

# 각 ligand 처리 시간 추정 (작은 DEG에서는 거의 0이었지만, 실제로는 계산됨)
# 보수적으로 ligand당 0.1초 가정 (작은 DEG 기준)
time_per_ligand_small <- 0.1  # 초
time_per_ligand_large <- time_per_ligand_small * scaling_factor * (top_targets_large / top_targets_small)
total_target_inference_time <- time_per_ligand_large * n_ligands

cat("3. Ligand-Target Inference (보수적 추정):\n")
cat("   - DEG 스케일링: sqrt(", round(deg_ratio, 1), ") = ", round(scaling_factor, 1), "배\n")
cat("   - Target 수 스케일링: ", top_targets_large, "/", top_targets_small, " = ", round(top_targets_large/top_targets_small, 2), "배\n")
cat("   - Ligand당 시간: ", round(time_per_ligand_large, 2), "초\n")
cat("   - 총 시간 (20 ligands): ", round(total_target_inference_time, 1), "초 (", round(total_target_inference_time/60, 1), "분)\n\n")

# 4. prepare_ligand_target_visualization 시간
# 이건 행렬 크기에 비례
viz_time_estimate <- (n_ligands * top_targets_large) / 1000  # 대략적인 추정
cat("4. Visualization Matrix Preparation:\n")
cat("   - 추정: ", round(viz_time_estimate, 1), "초\n\n")

# 5. Circos Plot (이미 최적화됨)
circos_time <- 0.1  # 초 (최적화 후)
cat("5. Circos Plot Generation:\n")
cat("   - 최적화 후: ", circos_time, "초\n\n")

# 총 시간 추정
total_time_seconds <- ligand_activity_time_large + total_target_inference_time + viz_time_estimate + circos_time
total_time_minutes <- total_time_seconds / 60

cat("=== 총 예상 시간 ===\n")
cat("보수적 추정: ", round(total_time_seconds, 0), "초 (", round(total_time_minutes, 1), "분)\n")
cat("코드 휴리스틱 기반: ", round((ligand_activity_time_large + heuristic_minutes*60 + viz_time_estimate + circos_time)/60, 1), "분\n\n")

cat("=== 참고 ===\n")
cat("- Ligand activity는 DEG 수에 비례하지만, 실제로는 ligand_target_matrix 크기가 더 중요\n")
cat("- Ligand-target inference는 sub-linear scaling (sqrt 가정)\n")
cat("- 실제 시간은 시스템 성능, 메모리 상태에 따라 달라질 수 있음\n")

