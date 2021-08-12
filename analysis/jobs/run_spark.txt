#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=08:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N spatial_proc
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_singleprocessing.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spatial_singleprocessing.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1.1_spatial_QC_reading/run_singleprocessing_spatial.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/visium_data/" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/objects/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/initial_qc/" \
        --force_filter;