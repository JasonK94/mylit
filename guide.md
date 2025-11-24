개요
scRNAseq / Spatial Transcriptomics 분석 툴 개발 토탈 패키지임.

각론
1. 테스트 환경:
cd /home/user3/GJC_KDW_250721 && Rscript 으로 스크립트를  실행하면, 패키지 의존성 문제는 없음. (renv에 모두 설치돼있음, start.R이 알아서 load해줌)
cursor가 직접 테스트 & 디버깅 루프를 반복하여 분석 결과를 /data/user3/sobj에 저장해두길 바람. overriding은 가능한 한 회피하여야 함.
2. 사용 데이터:
stroke PBMC full data: is5, 위치: /data/user3/sobj/IS6_sex_added_251110.qs
stroke PBMC downsampled data: is5s, 위치: /data/user3/sobj/IS6_sex_added_0.1x_251110.qs
mIBD GeoMx data: data_seurat, 위치: /data/user3/sobj/data_seurat_251104.qs
3. 변수 설명:
1) stroke
g3이 target variable("2", "1", "NA"), 환자는 $hos_no(8자리 숫자), 클러스터 정보는 $anno3.scvi(약 22가지)에 있다. $GEM, $set이 환자 외 random effect이지만, 환자는 23명, $GEM은 최대 12개, $set은 최대 3개이니 fixed effect 처리해도 무방함.
그러나 $hos_no<$GEM<$set 완전포함되어 동시에 넣을 시 공선성 문제가 발생하니 유의하고, 일부만 선정하여 사용할 것.
예시:
* gene_expression ~ g3 + anno3.scvi + (1|hos_no)
* (subset(anno3.scvi == "CD4+ T-cells") 등 subset 후에) gene_expression ~ g3 + (1|hos_no)
2) $emrid가 환자이며, 샘플수가 114개로 훨씬 적고, $drug(3가지 약제), $ck(조직 유형; PanCK positive or negative; TRUE, FALSE), $treatment ("pre", "post"), $response ("R", "NR")로 분석 고려사항은 많으며 데이터가 적어 더 까다로움.
그러나 양이 적으니 테스트하긴 빨라서 좋을 수도 있다. $emrid 별로 고유한 AOI(cell) 수가 2~8개 수준이니, 굳이 random effect로 고려 안 해도 될 수도 있음.
예시:
* simple: (drug, tissue type 별 subset 후에) gene_expression ~ treatment*response
* full model: gene_expression ~ treatment*ck*drug*response

4. 주의사항
1) 각 워크트리는, 각 브랜치에서만 작업내용을 커밋하고, 다른 워크트리, 특히 main이나 merger 브랜치에 영향을 주지 않도록 한다. _wt/에 각자 폴더로 존재함.
2) 각 브랜치 merge 시 docs, script가 꼬이지 않도록 mylit/docs/와 mylit/scripts에 하위 폴더로 따로 관리하여야 함.

5. 목적
1) 통계적으로 엄밀하면서도 강력한 분석(DEG 분석) 툴 개발: 주로 myR/analysis.R에 존재, analysis 브랜치, _wt/analysis 워크트리에서 개발중. 기존 main2 브랜치/워크트리였음.
2) DEG 분석 툴의 consensus를 구하는 툴 개발: 주로 /data/user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus에 소스코드, /data/user3/git_repo/_wt/deg-consensus/scripts/에 테스트 스크립트 존재, deg-consensus-dev브랜치, _wt/deg-consensus 워크트리에서 개발중.
3) 그룹을 잘 구분하는 signature 및 meta-learner 툴 개발: /data/user3/git_repo/_wt/fgs/myR/R/signature.R에서 소스코드 관리 및 개발 중, /data/user3/git_repo/_wt/fgs/myR/R/test.R에 일부 원치 않는 중복 함수 존재.
4) lds 브랜치/워크트리: limma-dream 후 SVA를 통해 공변량의 분산 설명력 측정 및 시각화. 251118-12:55 기준, limma-dream이 아닌 다른 deg list를 전달할 수 있어야 하지 않나 생각함. 아직 패치 착수 안함.
/data/user3/git_repo/_wt/lds/myR/R/lds_corrplot.R, /data/user3/git_repo/_wt/lds/myR/R/lds_08_heatmaps.R, /data/user3/git_repo/_wt/lds/myR/R/lds.R, /data/user3/git_repo/_wt/lds/scripts/lds
5) pt.umap 브랜치/워크트리: cell 단위인 scRNAseq를 환자 수준 데이터로 축약하여 dimension reduction을 통해 환자 그룹간 유사성, 차이를 분석하고자 함.
6) CCI(Nichenetr) / MILO / trajectory analysis 등 패키지 사용 편의성: milo, cci, pseudotime-dev 브랜치/워크트리
7) plotting 함수 편의성: plots-dev 브랜치/워크트리.

6. 추가 참조:
1) Docs
* 에이전트의 개발 & 디버깅 루틴을 위한 일반적인 내용은 /home/user3/data_user3/git_repo/mylit/docs/wt_test_general에 저장되어있음 
* 각 워크트리용 docs는 docs에 하위 폴더에, scripts는 scripts 하위 폴더에 정리.
2) scripts
* 각 하위 프로젝트의 script는 scripts/<하위브랜치명> 에 정리되어 있음.