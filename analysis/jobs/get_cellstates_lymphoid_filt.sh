#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=08:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N cell_states
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_filt_Lymphoid.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_filt_Lymphoid.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

mkdir ./results/ct_data_filt/Lymphoid;

$CONDA_PREFIX/bin/Rscript ./analysis/10_cellstates/get_filtered_states.R \
        --data_path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/Lymphoid/Lymphoid_states.rds" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.pdf" \
        --cell_class_excl "0,3,4" \
        --class_label "opt_state" \
        --start_res 0.1;

$CONDA_PREFIX/bin/Rscript ./analysis/utils/h5ad_mirrors.R \
        --scell_data "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.rds" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.h5ad";

$CONDA_PREFIX/bin/Rscript ./analysis/utils/sce_mirrors.R \
        --seurat_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.rds" \
        --sce_folder "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states_sce" \
        --assay "RNA" \
        --reduction "umap";

$CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/Lymphoid_states.rds" \
        --vars "orig.ident,opt_state" \
        --collapsed "yes" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data_filt/Lymphoid/ps_Lymphoid_states.rds" \
        --def_assay "RNA";
