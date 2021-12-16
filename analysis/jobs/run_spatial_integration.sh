#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=08:00:00
#PBS -l mem=200gb
#PBS -S /bin/bash
#PBS -N spatial_proc
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_integration.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_integration.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/integrate_objects.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/objects/" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides_qcint.pdf" \
        --def_assay "Spatial" \
        --batch_file "/beegfs/work/hd_wh241/MI_revisions/markers/visium_batch_ann.csv" \
        --default_resolution 0.5;

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

$CONDA_PREFIX/bin/Rscript ./analysis/utils/sce_mirrors.R \
        --seurat_file "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --sce_folder "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides_sce" \
        --assay "Spatial" \
        --reduction "umap_harmony";