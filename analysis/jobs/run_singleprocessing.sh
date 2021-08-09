#!/bin/bash

#SBATCH --job-name=processrnaseq
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/run_singleprocessing.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript run_singleprocessing.R \
        --folder \
        --path "/net/data.isilon/ag-saez/bq_shared/scellMI/raw_RNA/" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/processed_RNA/" \
        --out_fig_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/initial_qc/";