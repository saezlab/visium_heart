#!/bin/bash

#SBATCH --job-name=annstates
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/ann_states.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript map_states.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_wstates.rds" \
        --dictionary_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/state_annotations.txt" \
        --inherit_var "cell_type" \
        --ann_var "deconv_col";