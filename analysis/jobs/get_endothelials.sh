#!/bin/bash

#SBATCH --job-name=annotatestates
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_endothelials.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript get_cell_states.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/endothelial_states.rds" \
        --out_fig_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/endothelial_states.pdf" \
        --cell_class "endothelial_cells" \
        --class_label "cell_type";