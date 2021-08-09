#!/bin/bash

#SBATCH --job-name=cardfuncomics
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/plot_cardio_dea.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript plot_funcomics.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/cardiomyocytes_states.rds" \
        --dea_data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/cardio/cardiomyocytes_dea.rds" \
        --nfeats 10 \
        --ngenes_ORA 10 \
        --pvalue_ORA 0.15 \
        --pvalue 0.001 \
        --lfc 0.5 \
        --test_assays "RNA,progeny,dorothea" \
        --group_class "opt_state" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/cardio/cardiomyocytes";