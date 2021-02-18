#!/bin/bash

#SBATCH --job-name=annotate
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/annotate_integrateddata.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript annotate_object.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data.rds" \
        --dictionary_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/markers/test_annotations.txt" \
        --object_id "opt_clust_integrated" \
        --dictionary_id "opt_clust_integrated" \
        --new_variable "cell_type" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds";
$CONDA_PREFIX/bin/Rscript plot_knownmarkers.R \
        --id_label "cell_type" \
        --path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds" \
        --out_fig_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/";