좋습니다. 지금까지 정리한 **pathway 분석 파이프라인**을 실무 기준으로 바로 쓸 수 있게 표로 요약합니다.

---

# 🧭 Pathway 분석 전체 파이프라인 요약

| 단계               | 내용                          | 입력                  | 방법                      | 출력                 | 비고              |
| ---------------- | --------------------------- | ------------------- | ----------------------- | ------------------ | --------------- |
| Gene-level 분석    | DEG / mixed model / emmeans | expression + model  | edgeR / dream / emmeans | logFC, z, p        | gene-by-gene 통계 |
| Ranking 생성       | 전체 gene 순위화                 | z.ratio / t / logFC | 정렬                      | ranked vector      | GSEA용           |
| 유의 gene 선택       | cutoff 적용                   | padj, logFC         | threshold               | DEG list           | ORA용            |
| GSEA             | rank 기반 enrichment          | 전체 gene             | fgsea / GSEA            | NES, FDR           | **주력 방법**       |
| ORA              | 유의 gene enrichment          | DEG list            | hypergeometric          | OR, p              | 보조 방법           |
| DB 매핑            | gene→pathway                | gene set DB         | msigdbr 등               | pathway hit        | 여러 DB 병행        |
| Filtering        | significance 필터             | NES/FDR             | cutoff                  | hit pathway        | FDR 기준          |
| Redundancy 제거    | 중복 term 정리                  | pathway 결과          | simplify                | representative set | GO에서 중요         |
| Theme clustering | biological 묶음               | pathway             | manual/cluster          | theme              | 해석 단계           |
| 시각화              | plot                        | pathway 결과          | dotplot 등               | figure             | 논문용             |

---

# 🧪 ORA vs GSEA 비교

| 구분        | ORA            | GSEA        |
| --------- | -------------- | ----------- |
| 입력        | 유의 gene만       | **전체 gene** |
| 기반        | overlap        | ranking     |
| 통계        | hypergeometric | running sum |
| cutoff 영향 | 큼              | 없음          |
| 약한 신호     | 놓침             | **탐지 가능**   |
| omics 표준  | 보조             | **주력**      |
| 당신 분석     | 보조             | **권장**      |

---

# 🧬 주요 Pathway DB 비교

| DB           | 유형        | ORA | GSEA | 특징                | 권장도   |
| ------------ | --------- | --- | ---- | ----------------- | ----- |
| GO BP        | ontology  | ✔   | ✔    | 가장 넓음             | ★★★   |
| KEGG         | canonical | ✔   | ✔    | 경로 구조 명확          | ★★★   |
| Reactome     | curated   | ✔   | ✔    | signaling detail  | ★★★★  |
| Hallmark     | MSigDB    | —   | ✔    | 핵심 프로그램           | ★★★★★ |
| MSigDB C7    | immune    | —   | ✔    | immune signatures | ★★★★★ |
| MSigDB C6    | oncogenic | —   | ✔    | cancer signaling  | ★★★★  |
| WikiPathways | community | ✔   | ✔    | 보완용               | ★★    |

---

# 🧮 GSEA 파이프라인 요약

| 단계    | 내용                                |
| ----- | --------------------------------- |
| Input | gene + ranking score (z.ratio 권장) |
| 정렬    | statistic 기준 내림차순                 |
| 계산    | pathway별 running enrichment score |
| 정규화   | NES 계산                            |
| 검정    | permutation 기반 p                  |
| 보정    | FDR                               |
| 출력    | NES, p, FDR, leading edge         |
| 필터    | FDR < 0.25 (느슨) / <0.1 (엄격)       |

---

# 🧾 Ranking metric 선택 가이드

| 분석 방법       | 추천 ranking              |
| ----------- | ----------------------- |
| emmeans     | **z.ratio**             |
| mixed model | t-stat                  |
| edgeR       | signed logFC            |
| limma       | t                       |
| general     | sign(logFC) × -log10(p) |

---

# 📊 결과 시각화 세트

| Plot                 | 용도               |
| -------------------- | ---------------- |
| GSEA enrichment plot | pathway 대표 그림    |
| NES barplot          | drug 간 비교        |
| dotplot              | multi pathway 요약 |
| cnetplot             | gene–pathway 연결  |
| leading edge heatmap | core gene 패턴     |

---

원하면 다음 단계로
👉 **“drug별 emmeans 결과 → 자동 pathway 파이프라인 함수”**
👉 **IBD/immune/macrophage 전용 gene set 전략표**
👉 **논문 figure panel 템플릿**까지 바로 이어서 정리해줄게요.
