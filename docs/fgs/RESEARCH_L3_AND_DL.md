# Research Report: L3 Ensemble and Deep Learning in Disease Signature Discovery

## 1. L3 Ensemble (Stacking)

### 개념
L3 Ensemble(Stacking)은 여러 L2 모델(Base Learners)의 예측 결과(Probability 또는 Class)를 입력으로 받아, 최종 예측을 수행하는 메타 모델(Meta Learner)을 학습시키는 기법입니다.

### 가치와 효용성
*   **성능 향상**: 서로 다른 알고리즘(예: Tree-based vs Linear vs Neural Net)은 서로 다른 오류 패턴을 가집니다. 앙상블은 이러한 오류를 상쇄시켜 전반적인 예측 성능(AUC, Accuracy)을 향상시킬 수 있습니다.
*   **일반화**: 단일 모델보다 과적합(Overfitting)에 강건할 수 있습니다.

### 한계점 (L4, L5...?)
*   **해석력 저하**: 모델이 복잡해질수록 "어떤 유전자가 중요한가?"에 대한 답을 얻기 어려워집니다. 생물학적 마커 발굴에서는 해석력이 성능만큼이나 중요합니다.
*   **Diminishing Returns**: 레이어를 계속 쌓는다고 성능이 비례해서 오르지 않습니다. L3 정도가 실용적인 한계이며, 그 이상은 계산 비용 대비 효과가 미미합니다.

### 결론
TML 프레임워크에서 **L3 Ensemble은 충분히 가치가 있습니다**. `caretEnsemble` 등을 활용하여 구현할 수 있으며, 최종적인 "Diagnostic Score"를 산출하는 데 유용합니다. 하지만 그 이상의 레이어는 불필요합니다.

---

## 2. Deep Learning Approaches (DNN, VAE, Transformer)

### DNN (Deep Neural Networks)
*   **현황**: 이미 TML의 L2 모델 중 `mlp`, `mlpKerasDropout` 등이 DNN의 일종입니다.
*   **한계**: Tabular 데이터(유전자 발현량)에서는 Tree-based 모델(XGBoost, Ranger)이 DNN보다 성능이 좋거나 대등한 경우가 많습니다. 특히 샘플 수(환자 수)가 적은 scRNA-seq 코호트 분석에서는 DNN이 과적합되기 쉽습니다.

### VAE (Variational Autoencoders) - 예: scVI
*   **원리**: 유전자 발현 데이터를 저차원 잠재 공간(Latent Space)으로 압축하고 다시 복원합니다.
*   **활용**: 질병 상태와 정상 상태를 잠재 공간 상에서 분리하거나, 비선형적인 질병 궤적(Trajectory)을 모델링하는 데 탁월합니다.
*   **FGS/TML과의 비교**: scVI는 전체 유전자를 사용하므로 노이즈에 강하지만, 특정 "Signature Gene"을 딱 집어내기는 어렵습니다(Gene Ranking은 가능).

### Transformer (scGPT, Geneformer)
*   **최신 트렌드**: 대규모 scRNA-seq 데이터(수천만 개 세포)로 사전 학습(Pre-training)된 모델을 사용합니다.
*   **장점**: 적은 데이터로도 Fine-tuning을 통해 높은 성능을 낼 수 있습니다. 문맥(Context)을 이해하듯 유전자 간의 상호작용을 학습합니다.
*   **단점**: 모델이 매우 크고 무거우며, 학습 및 추론에 많은 자원이 필요합니다.

### 종합 비교 및 제언

| 접근법 | 해석력 (Interpretability) | 데이터 요구량 | 성능 (적은 데이터) | 성능 (많은 데이터) |
| :--- | :---: | :---: | :---: | :---: |
| **FGS/TML (Current)** | **High** (Explicit Genes) | Low ~ Medium | **High** | Medium |
| **DNN (Simple)** | Low | Medium | Low (Overfitting) | High |
| **VAE (scVI)** | Medium | Medium | Medium | High |
| **Transformer** | Low (Attention Map) | **High** (Pre-training) | **High** (Fine-tuning) | **Very High** |

**제언**:
1.  **L3 Ensemble 도입**: TML의 마지막 단계로 L3 Stacking을 추가하여 성능을 극대화하는 것은 좋은 전략입니다.
2.  **Cluster-specific Modeling**: 사용자가 제안한 대로, 세포 유형(Cluster)별로 별도의 모델을 만드는 것은 생물학적 맥락을 반영하는 가장 효과적인 방법입니다.
3.  **DL 도입**: 굳이 복잡한 Transformer를 도입하기보다는, 현재의 FGS/TML 프레임워크를 고도화(L3, Cluster-specific)하는 것이 "Disease Signature Discovery"라는 목적에 더 부합합니다.
