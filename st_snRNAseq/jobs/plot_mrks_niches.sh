#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N spatial_pseudobulk
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/plot_mrks_niches.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/plot_mrks_niches.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --id_label "opt_clust_integrated" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/" \
        --used_assay "Spatial";