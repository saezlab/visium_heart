#!/bin/bash

#SBATCH --job-name=nmf
#SBATCH -t 200:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/nmf_all.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
$CONDA_PREFIX/bin/Rscript integrated_nmf.R \
        --deconv_mats_folder "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/spotlight_mjr/" \
        --slide_files_folder "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_visium/" \
        --alias "slight" \
        --ks "3,4,5,6" \
        --out_slides_folder "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/spotlight_mjr_nmf/nmf_all_mjr_slight";