import anndata
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, mean_squared_error

def find_prognostic_genes_xgb(
    adata: anndata.AnnData,
    sample_label_col: str,
    prognosis_factor_col: str,
    cell_identity_col: str = None,
    target_cell_identity: str = None,
    use_pseudobulk: bool = True,
    n_top_genes: int = 50,
    test_size: float = 0.2,
    random_state: int = 42
) -> pd.DataFrame:
    """
    XGBoost를 사용하여 예후 인자와 관련된 유전자를 찾고 중요도 순으로 나열합니다.

    Args:
        adata (anndata.AnnData): 입력 AnnData 객체. adata.X는 유전자 발현 행렬이어야 합니다.
        sample_label_col (str): adata.obs에서 샘플/환자 식별자 컬럼 이름. 유사벌크 생성 시 사용됩니다.
        prognosis_factor_col (str): adata.obs에서 예후 인자 컬럼 이름.
                                     (수치형 또는 'high'/'low'와 같은 이진 범주형).
        cell_identity_col (str, optional): adata.obs에서 세포 식별 정보 컬럼 이름.
                                           특정 세포 유형만 분석 대상일 경우 사용. Defaults to None.
        target_cell_identity (str, optional): cell_identity_col에서 분석 대상이 되는 특정 세포 클러스터 레이블.
                                              cell_identity_col이 제공된 경우 필수. Defaults to None.
        use_pseudobulk (bool): 유사벌크 분석을 수행할지 여부.
                               True이면 sample_label_col을 기준으로 유사벌크를 생성합니다.
                               False이면 세포 단위로 분석합니다. Defaults to True.
        n_top_genes (int): 반환할 상위 유전자 개수. Defaults to 50.
        test_size (float): 테스트 데이터셋 비율. Defaults to 0.2.
        random_state (int): 재현성을 위한 난수 시드. Defaults to 42.

    Returns:
        pd.DataFrame: 'gene'와 'importance' 컬럼을 가진 데이터프레임.
                      예후 예측/추정에 중요한 유전자를 중요도 순으로 정렬하여 반환.
                      모델 성능(정확도 또는 MSE)도 함께 반환.
    """

    if cell_identity_col and target_cell_identity:
        print(f"Filtering for cell identity: {target_cell_identity}")
        adata_subset = adata[adata.obs[cell_identity_col] == target_cell_identity].copy()
        if adata_subset.n_obs == 0:
            raise ValueError(f"No cells found for identity: {target_cell_identity}")
    else:
        adata_subset = adata.copy()

    # 예후 인자 (Y) 준비
    # 환자 단위로 분석하려면, 환자별로 예후 인자가 하나여야 함.
    # 유사벌크 생성 시 환자별로 prognosis_factor_col 값이 유일한지 확인 필요.
    if use_pseudobulk:
        # 유사벌크 생성 (샘플별 유전자 발현 합계 또는 평균)
        # 여기서는 합계를 사용. 필요시 평균 등으로 변경 가능.
        # 또한, 각 샘플에 대한 prognosis_factor_col 값이 유일해야 함.
        # sample_label_col을 기준으로 그룹화하고, 첫 번째 prognosis_factor_col 값을 사용한다고 가정.
        # 실제 데이터에서는 이 부분을 더 견고하게 만들어야 함 (예: 환자 메타데이터와 병합)

        if sample_label_col not in adata_subset.obs.columns:
            raise ValueError(f"Sample label column '{sample_label_col}' not found in adata.obs.")
        if prognosis_factor_col not in adata_subset.obs.columns:
            raise ValueError(f"Prognosis factor column '{prognosis_factor_col}' not found in adata.obs.")

        # 각 샘플에 대한 예후 인자 추출 (중복 제거)
        prognosis_data = adata_subset.obs[[sample_label_col, prognosis_factor_col]].drop_duplicates(subset=[sample_label_col])
        prognosis_data = prognosis_data.set_index(sample_label_col)

        # 유사벌크 생성
        pseudobulk_list = []
        for sample_id in adata_subset.obs[sample_label_col].unique():
            sample_cells = adata_subset[adata_subset.obs[sample_label_col] == sample_id, :]
            if sample_cells.n_obs > 0:
                # scipy.sparse.csr_matrix인 경우 .sum(axis=0) 후 .A1으로 변환
                if hasattr(sample_cells.X, "sum"):
                    pseudobulk_expression = sample_cells.X.sum(axis=0)
                    if isinstance(pseudobulk_expression, np.matrix): # scipy.sparse.matrix의 sum 결과
                         pseudobulk_expression = pseudobulk_expression.A1
                    elif not isinstance(pseudobulk_expression, np.ndarray): # 다른 sparse 타입의 sum 결과
                         pseudobulk_expression = pseudobulk_expression.toarray().ravel()
                else: # numpy.ndarray인 경우
                    pseudobulk_expression = sample_cells.X.sum(axis=0)

                df_temp = pd.DataFrame([pseudobulk_expression], columns=adata_subset.var_names, index=[sample_id])
                pseudobulk_list.append(df_temp)

        if not pseudobulk_list:
            raise ValueError("Pseudobulk creation resulted in an empty list. Check sample labels and data.")

        X_df = pd.concat(pseudobulk_list)
        # X_df의 인덱스와 prognosis_data의 인덱스를 기준으로 정렬 및 필터링
        common_samples = X_df.index.intersection(prognosis_data.index)
        if common_samples.empty:
            raise ValueError("No common samples found between pseudobulk data and prognosis data. Check sample IDs.")

        X_df = X_df.loc[common_samples]
        y_series = prognosis_data.loc[common_samples, prognosis_factor_col]
        X = X_df.values
        y = y_series.values

    else: # 세포 단위 분석
        X = adata_subset.X.toarray() if hasattr(adata_subset.X, "toarray") else adata_subset.X # 밀집 행렬로 변환
        y = adata_subset.obs[prognosis_factor_col].values

    # Y 데이터 타입 확인 및 변환
    task_type = None
    if pd.api.types.is_numeric_dtype(y):
        # 값이 두 개 뿐인 수치형 (예: 0, 1)은 분류로 처리 가능
        if len(np.unique(y)) == 2:
            print("Prognosis factor is numeric with 2 unique values. Treating as classification.")
            task_type = "classification"
            le = LabelEncoder()
            y = le.fit_transform(y)
        elif len(np.unique(y)) < 2:
             raise ValueError("Prognosis factor has less than 2 unique values. Cannot perform analysis.")
        else:
            print("Prognosis factor is numeric. Treating as regression.")
            task_type = "regression"
    elif pd.api.types.is_string_dtype(y) or pd.api.types.is_categorical_dtype(y):
        if len(np.unique(y)) == 2:
            print("Prognosis factor is categorical with 2 unique values. Treating as binary classification.")
            task_type = "classification"
            le = LabelEncoder()
            y = le.fit_transform(y) # 예: 'high'/'low' -> 0/1
        elif len(np.unique(y)) < 2:
             raise ValueError("Prognosis factor has less than 2 unique values. Cannot perform analysis.")
        else:
            # 다중 클래스 분류는 여기서는 지원하지 않음 (간단하게 이진 분류 또는 회귀에 초점)
            raise ValueError("Multi-class classification is not directly supported by this simplified function. "
                             "Please ensure prognosis factor is numeric or binary categorical (e.g., 'high'/'low').")
    else:
        raise ValueError("Prognosis factor data type is not recognized or suitable for XGBoost.")

    if X.shape[0] <= 1 or len(y) <= 1 :
        raise ValueError("Not enough samples to perform train/test split and XGBoost training.")
    if X.shape[0] != len(y):
        raise ValueError(f"Shape mismatch: X has {X.shape[0]} samples, y has {len(y)} samples.")


    # 데이터 분할
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state, stratify=y if task_type == "classification" else None)

    # XGBoost 모델 설정 및 학습
    model_performance = None
    if task_type == "classification":
        model = xgb.XGBClassifier(objective='binary:logistic', eval_metric='logloss', use_label_encoder=False, random_state=random_state)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        model_performance = {"metric": "accuracy", "value": accuracy}
        print(f"XGBoost Classifier Accuracy: {accuracy:.4f}")
    elif task_type == "regression":
        model = xgb.XGBRegressor(objective='reg:squarederror', eval_metric='rmse', random_state=random_state)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mse = mean_squared_error(y_test, y_pred)
        model_performance = {"metric": "mse", "value": mse}
        print(f"XGBoost Regressor MSE: {mse:.4f}")
    else:
        raise ValueError("Undefined task type.") # Should not happen

    # 특성 중요도 추출
    feature_importances = model.feature_importances_
    gene_names = adata_subset.var_names if not use_pseudobulk else X_df.columns

    importance_df = pd.DataFrame({
        'gene': gene_names,
        'importance': feature_importances
    }).sort_values(by='importance', ascending=False)

    return importance_df.head(n_top_genes), model_performance

# --- 사용 예시 (실제 실행을 위해서는 AnnData 객체와 데이터가 필요합니다) ---
# import scanpy as sc
#
# # 1. AnnData 객체 로드 또는 생성 (예시)
# # adata = sc.read_h5ad("path_to_your_data.h5ad")
# # 가상 데이터 생성
# n_cells = 200
# n_genes = 1000
# n_samples = 10
#
# X_data = np.random.rand(n_cells, n_genes)
# obs_data = pd.DataFrame({
#     'sample_id': [f'sample_{i%n_samples + 1}' for i in range(n_cells)],
#     'cell_type': np.random.choice(['A', 'B', 'C'], size=n_cells),
#     # 각 샘플별로 예후 인자 할당 (이진)
#     'prognosis_binary': ['high' if (i%n_samples + 1) % 2 == 0 else 'low' for i in range(n_cells)],
#     # 각 샘플별로 예후 인자 할당 (수치)
#     'prognosis_numeric': [(i%n_samples + 1) * 0.5 + np.random.normal(0,0.1) for i in range(n_cells)],
# })
# var_data = pd.DataFrame(index=[f'gene_{j}' for j in range(n_genes)])
#
# adata_example = anndata.AnnData(X=X_data, obs=obs_data, var=var_data)
# adata_example.obs['prognosis_binary_per_sample'] = adata_example.obs.groupby('sample_id')['prognosis_binary'].transform('first')
# adata_example.obs['prognosis_numeric_per_sample'] = adata_example.obs.groupby('sample_id')['prognosis_numeric'].transform('first')
#
# # 예시 1: 유사벌크, 이진 예후 인자, 특정 세포 타입 필터링
# try:
#     print("\n--- Example 1: Pseudobulk, Binary Prognosis, Cell Type A ---")
#     top_genes_pb_binary, perf_pb_binary = find_prognostic_genes_xgb(
#         adata=adata_example,
#         sample_label_col='sample_id',
#         prognosis_factor_col='prognosis_binary_per_sample', # 환자별로 값이 동일해야 함
#         cell_identity_col='cell_type',
#         target_cell_identity='A',
#         use_pseudobulk=True,
#         n_top_genes=10
#     )
#     print("Top prognostic genes (Pseudobulk, Binary):")
#     print(top_genes_pb_binary)
#     print("Model Performance:", perf_pb_binary)
# except Exception as e:
#     print(f"Error in Example 1: {e}")
#
# # 예시 2: 세포 단위, 수치 예후 인자
# try:
#     print("\n--- Example 2: Cell-level, Numeric Prognosis ---")
#     # 세포 단위 분석을 위해서는 각 세포마다 예후 인자 값이 있어야 합니다.
#     # 여기서는 임의로 'prognosis_numeric'을 사용 (실제로는 의미있는 세포 단위 예후 인자여야 함)
#     adata_example.obs['prognosis_numeric_cell_level'] = np.random.rand(n_cells) * 10
#     top_genes_cell_numeric, perf_cell_numeric = find_prognostic_genes_xgb(
#         adata=adata_example,
#         sample_label_col='sample_id', # 세포 단위 분석 시에는 덜 중요하지만, 구조상 필요
#         prognosis_factor_col='prognosis_numeric_cell_level',
#         use_pseudobulk=False,
#         n_top_genes=10
#     )
#     print("Top prognostic genes (Cell-level, Numeric):")
#     print(top_genes_cell_numeric)
#     print("Model Performance:", perf_cell_numeric)
# except Exception as e:
#     print(f"Error in Example 2: {e}")
#
# # 예시 3: 유사벌크, 수치 예후 인자
# try:
#     print("\n--- Example 3: Pseudobulk, Numeric Prognosis ---")
#     top_genes_pb_numeric, perf_pb_numeric = find_prognostic_genes_xgb(
#         adata=adata_example,
#         sample_label_col='sample_id',
#         prognosis_factor_col='prognosis_numeric_per_sample', # 환자별로 값이 동일해야 함
#         use_pseudobulk=True,
#         n_top_genes=10
#     )
#     print("Top prognostic genes (Pseudobulk, Numeric):")
#     print(top_genes_pb_numeric)
#     print("Model Performance:", perf_pb_numeric)
# except Exception as e:
#     print(f"Error in Example 3: {e}")



import anndata
import numpy as np
import pandas as pd
from collections import defaultdict
# find_prognostic_genes_xgb 함수가 정의되어 있다고 가정합니다.
# 만약 다른 파일에 있다면 from your_module import find_prognostic_genes_xgb 와 같이 가져옵니다.

def find_robust_prognostic_genes_with_downsampling(
    adata: anndata.AnnData,
    sample_label_col: str,
    prognosis_factor_col: str,
    n_iterations: int = 10,
    downsampling_fraction: float = 0.8, # 각 반복에서 사용할 데이터의 비율
    use_pseudobulk: bool = True,
    cell_identity_col: str = None,
    target_cell_identity: str = None,
    n_final_top_genes: int = 50,
    random_state_start: int = 0 # 각 반복에 다른 random_state를 주기 위한 시작점
) -> pd.DataFrame:
    """
    다운샘플링과 반복적인 XGBoost 실행을 통해 안정적으로 예후 관련 유전자를 찾습니다.

    Args:
        adata (anndata.AnnData): 입력 AnnData 객체.
        sample_label_col (str): adata.obs에서 샘플/환자 식별자 컬럼 이름.
        prognosis_factor_col (str): adata.obs에서 예후 인자 컬럼 이름.
        n_iterations (int): 다운샘플링 및 XGBoost 실행 반복 횟수. Defaults to 10.
        downsampling_fraction (float): 각 반복에서 사용할 데이터(세포 또는 유사벌크 샘플)의 비율.
                                      Defaults to 0.8 (80%).
        use_pseudobulk (bool): 유사벌크 분석을 수행할지 여부. Defaults to True.
        cell_identity_col (str, optional): 세포 식별 정보 컬럼. Defaults to None.
        target_cell_identity (str, optional): 분석 대상 특정 세포 클러스터 레이블. Defaults to None.
        n_final_top_genes (int): 최종적으로 반환할 상위 유전자 개수. Defaults to 50.
        random_state_start (int): 반복마다 다른 random_state를 생성하기 위한 시작 시드. Defaults to 0.

    Returns:
        pd.DataFrame: 'gene', 'average_importance', 'frequency' 컬럼을 가진 데이터프레임.
                      여러 반복에 걸쳐 집계된 유전자 중요도 및 선택 빈도 순으로 정렬.
    """

    all_gene_importances = defaultdict(list)
    gene_selection_counts = defaultdict(int)
    all_genes_set = set(adata.var_names) # 모든 유전자 이름

    # 초기 필터링 (한 번만 수행)
    if cell_identity_col and target_cell_identity:
        print(f"Initial filtering for cell identity: {target_cell_identity}")
        adata_filtered = adata[adata.obs[cell_identity_col] == target_cell_identity].copy()
        if adata_filtered.n_obs == 0:
            raise ValueError(f"No cells found for identity: {target_cell_identity}")
    else:
        adata_filtered = adata.copy()

    if use_pseudobulk:
        # 유사벌크를 미리 만들고, 유사벌크 샘플을 다운샘플링
        # (주의: 이 부분은 find_prognostic_genes_xgb 내부 로직과 유사하게 처리하거나,
        #  find_prognostic_genes_xgb가 AnnData 대신 DataFrame/ndarray를 받도록 수정 필요)

        # 간단화된 유사벌크 생성 (find_prognostic_genes_xgb의 로직 일부 활용)
        prognosis_data_full = adata_filtered.obs[[sample_label_col, prognosis_factor_col]].drop_duplicates(subset=[sample_label_col])
        prognosis_data_full = prognosis_data_full.set_index(sample_label_col)

        pseudobulk_list_full = []
        for sample_id_full in adata_filtered.obs[sample_label_col].unique():
            sample_cells_full = adata_filtered[adata_filtered.obs[sample_label_col] == sample_id_full, :]
            if sample_cells_full.n_obs > 0:
                if hasattr(sample_cells_full.X, "sum"):
                    pseudobulk_expression_full = sample_cells_full.X.sum(axis=0)
                    if isinstance(pseudobulk_expression_full, np.matrix):
                         pseudobulk_expression_full = pseudobulk_expression_full.A1
                    elif not isinstance(pseudobulk_expression_full, np.ndarray):
                         pseudobulk_expression_full = pseudobulk_expression_full.toarray().ravel()
                else:
                    pseudobulk_expression_full = sample_cells_full.X.sum(axis=0)
                df_temp_full = pd.DataFrame([pseudobulk_expression_full], columns=adata_filtered.var_names, index=[sample_id_full])
                pseudobulk_list_full.append(df_temp_full)

        if not pseudobulk_list_full:
            raise ValueError("Initial pseudobulk creation resulted in an empty list.")
        X_df_full = pd.concat(pseudobulk_list_full)
        # X_df_full의 인덱스와 prognosis_data_full의 인덱스를 기준으로 정렬 및 필터링
        common_samples_full = X_df_full.index.intersection(prognosis_data_full.index)
        if common_samples_full.empty:
            raise ValueError("No common samples found between full pseudobulk data and prognosis data.")

        X_df_full = X_df_full.loc[common_samples_full]
        y_series_full = prognosis_data_full.loc[common_samples_full, prognosis_factor_col]

        if X_df_full.shape[0] < 2: # 다운샘플링 및 train/test split을 하기에 너무 적은 샘플
             raise ValueError(f"Not enough pseudobulk samples ({X_df_full.shape[0]}) for downsampling and analysis.")
        
        original_indices = X_df_full.index

    else: # 세포 단위 분석
        if adata_filtered.n_obs < 2: # 다운샘플링 및 train/test split을 하기에 너무 적은 세포
             raise ValueError(f"Not enough cells ({adata_filtered.n_obs}) for downsampling and analysis.")
        original_indices = adata_filtered.obs.index


    for i in range(n_iterations):
        current_random_state = random_state_start + i
        print(f"\nIteration {i+1}/{n_iterations} with random_state {current_random_state}")

        # 데이터 다운샘플링
        n_samples_to_select = int(len(original_indices) * downsampling_fraction)
        if n_samples_to_select < 2 : # 최소 2개는 있어야 train/test split 가능
            print(f"Skipping iteration {i+1} due to insufficient samples after downsampling ({n_samples_to_select}).")
            continue

        np.random.seed(current_random_state)
        downsampled_indices = np.random.choice(original_indices, size=n_samples_to_select, replace=False)

        if use_pseudobulk:
            X_downsampled_df = X_df_full.loc[downsampled_indices]
            y_downsampled_series = y_series_full.loc[downsampled_indices]

            # 유사벌크 데이터를 find_prognostic_genes_xgb가 처리할 수 있는 형태로 변환
            # find_prognostic_genes_xgb가 DataFrame과 Series를 직접 받을 수 있도록 수정하거나,
            # 여기서 임시 AnnData 객체를 만들어 전달해야 함.
            # 여기서는 임시 AnnData를 만드는 방식을 가정.
            # (더 효율적인 방법은 find_prognostic_genes_xgb를 수정하는 것)
            temp_obs = pd.DataFrame({prognosis_factor_col: y_downsampled_series,
                                     sample_label_col: y_downsampled_series.index # 임시로 sample_label_col에 인덱스 사용
                                     }, index=y_downsampled_series.index)

            # X_downsampled_df.columns (유전자 이름)을 adata_filtered.var_names 와 일치시켜야 함.
            # X_df_full 생성 시 이미 var_names를 사용했으므로 일치함.
            temp_adata = anndata.AnnData(X=X_downsampled_df.values,
                                         obs=temp_obs,
                                         var=pd.DataFrame(index=X_downsampled_df.columns))
            # find_prognostic_genes_xgb는 sample_label_col을 obs에서 찾으려 함.
            # use_pseudobulk=True로 호출하지만, 이미 데이터가 유사벌크이므로,
            # find_prognostic_genes_xgb 내부의 유사벌크 로직은 사실상 현재 데이터를 그대로 사용하게 됨.
            # (이 부분이 다소 복잡하므로, find_prognostic_genes_xgb를 좀 더 유연하게 만드는 것이 좋음)
            # 여기서는 find_prognostic_genes_xgb가 이미 pseudobulk된 데이터를 받을 때,
            # sample_label_col을 통해 y값을 매칭한다고 가정하고 진행.
            # 이 경우 find_prognostic_genes_xgb의 유사벌크 생성부분은 건너뛰도록 수정이 필요할 수 있음
            # 또는, 가장 간단하게는 use_pseudobulk=False로 하고, X, y를 직접 전달하는 형태로 함수를 수정.

            # find_prognostic_genes_xgb를 직접 호출하는 대신, 그 핵심 로직을 여기에 통합하거나,
            # find_prognostic_genes_xgb가 다운샘플링된 DataFrame을 직접 받을 수 있게 수정.
            # 아래는 find_prognostic_genes_xgb를 그대로 사용한다고 가정하고, 임시 AnnData를 사용.
            # 단, 이 경우 find_prognostic_genes_xgb의 use_pseudobulk=True는
            # 내부에서 다시 유사벌크를 만들려고 시도할 수 있으므로 주의.
            # 임시 해결: find_prognostic_genes_xgb에 X, y를 직접 전달하는 옵션 추가 또는
            # 여기서 직접 XGBoost 학습 로직을 실행. 여기서는 후자를 택해 간략화.

            X_iter = X_downsampled_df.values
            y_iter = y_downsampled_series.values
            iter_gene_names = X_downsampled_df.columns

        else: # 세포 단위 다운샘플링
            adata_downsampled = adata_filtered[downsampled_indices, :].copy()
            if adata_downsampled.n_obs == 0:
                print(f"Skipping iteration {i+1} due to empty data after cell downsampling.")
                continue
            # find_prognostic_genes_xgb 함수를 직접 사용.
            # 이 경우 sample_label_col은 find_prognostic_genes_xgb 내부에서만 사용됨.
            try:
                # 이 부분은 find_prognostic_genes_xgb의 핵심 로직을 가져와 직접 실행하는 것이 더 깔끔함
                # 지금은 함수 호출로 대체 (내부에서 use_pseudobulk=False로 작동)
                top_genes_df_iter, _ = find_prognostic_genes_xgb(
                    adata=adata_downsampled, # 다운샘플링된 AnnData
                    sample_label_col=sample_label_col, # 그대로 전달
                    prognosis_factor_col=prognosis_factor_col,
                    use_pseudobulk=False, # 세포 단위로 이미 준비됨
                    n_top_genes=adata_downsampled.n_vars, # 모든 유전자 중요도 받기
                    random_state=current_random_state
                )
                # 모든 유전자에 대해 중요도 기록
                for _, row in top_genes_df_iter.iterrows():
                    all_gene_importances[row['gene']].append(row['importance'])
                    if row['importance'] > 0: # 중요도가 0보다 큰 경우에만 선택된 것으로 간주
                        gene_selection_counts[row['gene']] += 1

            except ValueError as e:
                print(f"Skipping iteration {i+1} due to error in find_prognostic_genes_xgb: {e}")
                continue
            except Exception as e:
                print(f"An unexpected error occurred in iteration {i+1}: {e}")
                continue
            continue # 세포 단위 분석 후 다음 반복으로


        # --- 유사벌크 데이터에 대해 직접 XGBoost 로직 실행 (find_prognostic_genes_xgb 대신) ---
        # (find_prognostic_genes_xgb를 재사용하려면 입력 형식 맞추는 부분이 복잡해져서 직접 구현)
        if use_pseudobulk: # 이 블록은 use_pseudobulk=True일 때만 실행
            task_type_iter = None
            y_transformed_iter = y_iter.copy() # 원본 y_iter를 변경하지 않기 위해 복사

            if pd.api.types.is_numeric_dtype(y_transformed_iter):
                if len(np.unique(y_transformed_iter)) == 2:
                    task_type_iter = "classification"
                    le_iter = LabelEncoder()
                    y_transformed_iter = le_iter.fit_transform(y_transformed_iter)
                elif len(np.unique(y_transformed_iter)) < 2:
                    print(f"Skipping iteration {i+1} due to <2 unique prognosis values in subset.")
                    continue
                else:
                    task_type_iter = "regression"
            elif pd.api.types.is_string_dtype(y_transformed_iter) or pd.api.types.is_categorical_dtype(y_transformed_iter):
                if len(np.unique(y_transformed_iter)) == 2:
                    task_type_iter = "classification"
                    le_iter = LabelEncoder()
                    y_transformed_iter = le_iter.fit_transform(y_transformed_iter)
                elif len(np.unique(y_transformed_iter)) < 2:
                    print(f"Skipping iteration {i+1} due to <2 unique prognosis values in subset.")
                    continue
                else:
                    print(f"Skipping iteration {i+1} due to multi-class prognosis factor in subset.")
                    continue # 다중 클래스 건너뛰기
            else:
                print(f"Skipping iteration {i+1} due to unrecognized prognosis factor type in subset.")
                continue

            if X_iter.shape[0] < 2 or len(y_transformed_iter) < 2: # train/test split에 충분한 샘플 확인
                print(f"Skipping iteration {i+1} due to insufficient samples for train/test split after downsampling.")
                continue


            X_train_iter, _, y_train_iter, _ = train_test_split(
                X_iter, y_transformed_iter, test_size=0.2, random_state=current_random_state, # test_size는 find_prognostic_genes_xgb와 동일하게
                stratify=y_transformed_iter if task_type_iter == "classification" else None
            )

            if X_train_iter.shape[0] == 0:
                print(f"Skipping iteration {i+1} due to empty training set after split.")
                continue

            model_iter = None
            if task_type_iter == "classification":
                model_iter = xgb.XGBClassifier(objective='binary:logistic', eval_metric='logloss', use_label_encoder=False, random_state=current_random_state)
            elif task_type_iter == "regression":
                model_iter = xgb.XGBRegressor(objective='reg:squarederror', eval_metric='rmse', random_state=current_random_state)

            if model_iter is None:
                print(f"Skipping iteration {i+1} as model type could not be determined.")
                continue

            try:
                model_iter.fit(X_train_iter, y_train_iter) # 테스트셋 없이 전체 다운샘플링 데이터로 학습 후 중요도만 추출
                                                       # 또는, train/test split 후 학습 데이터로만 학습
            except Exception as e:
                print(f"Error fitting XGBoost model in iteration {i+1}: {e}")
                continue

            importances_iter = model_iter.feature_importances_
            for gene_idx, importance_val in enumerate(importances_iter):
                gene_name = iter_gene_names[gene_idx]
                all_gene_importances[gene_name].append(importance_val)
                if importance_val > 0:
                    gene_selection_counts[gene_name] += 1
        # --- 유사벌크 XGBoost 로직 끝 ---


    # 모든 반복으로부터 결과 집계
    aggregated_results = []
    for gene in all_genes_set: # adata.var_names 대신 all_genes_set 사용
        importances = all_gene_importances.get(gene, [0]) # 한 번도 선택 안된 유전자는 중요도 0으로 처리
        avg_importance = np.mean(importances) if importances else 0
        frequency = gene_selection_counts.get(gene, 0) / n_iterations # 선택된 빈도
        aggregated_results.append({
            'gene': gene,
            'average_importance': avg_importance,
            'frequency': frequency,
            'sum_importance': np.sum(importances) # 또는 합계 중요도
        })

    results_df = pd.DataFrame(aggregated_results)
    # 평균 중요도와 빈도를 기준으로 정렬 (필요에 따라 다른 기준으로 변경 가능)
    results_df = results_df.sort_values(by=['average_importance', 'frequency'], ascending=[False, False])

    return results_df.head(n_final_top_genes)


# --- 사용 예시 (실제 실행을 위해서는 AnnData 객체와 데이터, find_prognostic_genes_xgb 함수가 필요합니다) ---
# import scanpy as sc
#
# # 가상 데이터 생성 (이전과 동일)
# n_cells = 300
# n_genes = 500
# n_samples = 15
#
# X_data = np.random.rand(n_cells, n_genes) + np.random.normal(0,0.2, (n_cells,n_genes)) # 약간의 노이즈 추가
# X_data = np.abs(X_data) # 발현량은 양수
#
# # 일부 유전자에 신호 추가
# prognostic_genes_indices = np.random.choice(n_genes, 10, replace=False) # 10개의 예후 관련 유전자
#
# obs_list = []
# for i in range(n_cells):
#     sample_num = (i % n_samples) + 1
#     is_high_prog = (sample_num % 3 == 0) # 샘플 3, 6, 9... 가 'high'
#     prog_binary = 'high' if is_high_prog else 'low'
#     prog_numeric = 10 + sample_num * 0.2 + np.random.normal(0,0.1) if is_high_prog else 5 - sample_num*0.1 + np.random.normal(0,0.1)
#
#     obs_list.append({
#         'sample_id': f'sample_{sample_num}',
#         'cell_type': np.random.choice(['TypeX', 'TypeY']),
#         'prognosis_binary_per_sample': prog_binary,
#         'prognosis_numeric_per_sample': prog_numeric,
#     })
#     # 신호가 있는 유전자 발현량 조작
#     if is_high_prog:
#         X_data[i, prognostic_genes_indices] += np.random.uniform(0.5, 1.5, len(prognostic_genes_indices))
#     else:
#         X_data[i, prognostic_genes_indices] -= np.random.uniform(0.1, 0.5, len(prognostic_genes_indices))
# X_data = np.maximum(X_data, 0) # 발현량이 음수가 되지 않도록
#
# obs_data = pd.DataFrame(obs_list)
# var_data = pd.DataFrame(index=[f'gene_{j}' for j in range(n_genes)])
#
# adata_robust_example = anndata.AnnData(X=X_data, obs=obs_data, var=var_data)
#
# print(f"True prognostic gene indices: {prognostic_genes_indices}")
# print(f"True prognostic gene names: {[f'gene_{j}' for j in prognostic_genes_indices]}")
#
# # 예시 1: 유사벌크, 이진 예후, 5회 반복
# try:
#     print("\n--- Robust Example 1: Pseudobulk, Binary Prognosis, 5 Iterations ---")
#     robust_genes_pb_binary = find_robust_prognostic_genes_with_downsampling(
#         adata=adata_robust_example,
#         sample_label_col='sample_id',
#         prognosis_factor_col='prognosis_binary_per_sample',
#         n_iterations=5, # 반복 횟수 줄여서 테스트 시간 단축
#         downsampling_fraction=0.7,
#         use_pseudobulk=True,
#         # cell_identity_col='cell_type', # 필요시 특정 세포 타입만 사용
#         # target_cell_identity='TypeX',
#         n_final_top_genes=20,
#         random_state_start=100
#     )
#     print("Robust Top Prognostic Genes (Pseudobulk, Binary):")
#     print(robust_genes_pb_binary)
# except Exception as e:
#     print(f"Error in Robust Example 1: {e}")
#     import traceback
#     traceback.print_exc()
#
#
# # 예시 2: 세포 단위, 수치 예후, 3회 반복 (세포 단위는 시간이 오래 걸릴 수 있음)
# # 세포 단위 분석을 위해 각 세포에 고유한 예후 인자 필요 (여기서는 샘플 레벨 값을 그대로 사용)
# adata_robust_example.obs['prognosis_numeric_cell_level'] = adata_robust_example.obs['prognosis_numeric_per_sample']
# try:
#     print("\n--- Robust Example 2: Cell-level, Numeric Prognosis, 3 Iterations ---")
#     robust_genes_cell_numeric = find_robust_prognostic_genes_with_downsampling(
#         adata=adata_robust_example,
#         sample_label_col='sample_id', # 세포 단위에서는 다운샘플링 인덱싱에만 활용
#         prognosis_factor_col='prognosis_numeric_cell_level',
#         n_iterations=3, # 반복 횟수 줄임
#         downsampling_fraction=0.1, # 세포가 많으므로 샘플링 비율을 낮춤
#         use_pseudobulk=False,
#         n_final_top_genes=20,
#         random_state_start=200
#     )
#     print("Robust Top Prognostic Genes (Cell-level, Numeric):")
#     print(robust_genes_cell_numeric)
# except Exception as e:
#     print(f"Error in Robust Example 2: {e}")
#     import traceback
#     traceback.print_exc()