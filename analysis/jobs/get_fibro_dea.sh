#!/bin/bash

#SBATCH --job-name=cardfuncomics
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_fibro_dea.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript estimate_dea.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/fibroblasts_states.rds" \
        --out_df "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/fibro/fibroblasts_dea.rds" \
        --test_assays "RNA,progeny,dorothea,ecm_scores" \
        --group_class "opt_state" \
        --lfc "0.05" \
        --only_pos "yes";