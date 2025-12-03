runid=run10


Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe_validate.R \
  --config config/config_complete2.csv \
  --run_id $runid

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe1_read_demulti.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 0 \
  --output_step 1

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe2_nmz_clustering.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 1 \
  --output_step 2

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe3_ambient_removal.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 2 \
  --output_step 3

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe4_sctransform.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 3 \
  --output_step 4

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe5_doubletfinder.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 4 \
  --output_step 5

Rscript /data/user3/git_repo/_wt/pipe/scripts/pipe6_integration.R \
  --config config/config_complete2.csv \
  --run_id $runid \
  --input_step 5 \
  --output_step 6 \
  --method scVI


