#!/bin/bash

#SBATCH --job-name=cardfuncomics
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_cardio_funcomics.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript add_funcomics.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/cardiomyocytes_states.rds" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/cardiomyocytes_states.rds" \
        --group_class "opt_state" \
        --gene_set_collection "not";