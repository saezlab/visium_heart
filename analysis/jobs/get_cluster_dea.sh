#!/bin/bash

#SBATCH --job-name=getclusterdegenes
#SBATCH -t 300:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_cluster_dea.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript estimate_dea.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data.rds" \
        --out_df "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_degs.rds" \
        --test_assays "RNA" \
        --group_class "opt_clust_integrated" \
        --lfc "0.05" \
        --only_pos "yes";