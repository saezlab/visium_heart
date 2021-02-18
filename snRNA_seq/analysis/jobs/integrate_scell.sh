#!/bin/bash

#SBATCH --job-name=integrate
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/integrate_scell.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript integrate_objects.R \
        --path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/processed_RNA/" \
        --out_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data.rds" \
        --out_fig_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_qcint.pdf";
$CONDA_PREFIX/bin/Rscript plot_knownmarkers.R \
        --id_label "opt_clust_integrated" \
        --path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data.rds" \
        --out_fig_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/";