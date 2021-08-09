#!/bin/bash

#SBATCH --job-name=deconv
#SBATCH -t 300:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/deconv_slight.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript run_spotlight.R \
        --scell_data "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds" \
        --ident_var "cell_type" \
        --folder \
        --data_path "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_visium/" \
        --out_fig_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/spotlight_mjr/" \
        --out_file "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/spotlight_mjr/";