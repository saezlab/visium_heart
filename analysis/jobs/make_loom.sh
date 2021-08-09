#!/bin/bash

#SBATCH --job-name=looming
#SBATCH -t 100:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/make_loom.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript ./utils/loom_mirrors.R \
        --scell_data "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_wstates.rds" \
        --out_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_wstates.loom";