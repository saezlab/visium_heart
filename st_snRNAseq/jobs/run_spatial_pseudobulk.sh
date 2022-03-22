#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N spatial_pseudobulk
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_pseudobulk.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_pseudobulk.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --vars "orig.ident" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/ps_integrated_slides.rds" \
        --def_assay "Spatial";

$CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --vars "orig.ident,opt_clust_integrated" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/ps_integrated_slides_niches.rds" \
        --collapsed "yes" \
        --def_assay "Spatial";