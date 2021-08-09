#!/bin/bash

#SBATCH --job-name=getpseudo
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_pseudobulk_atlas.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript ./get_pseudobulk_profiles.R \
        --path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds" \
        --vars "cell_type" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/ps_integrated_data_annotated.rds";