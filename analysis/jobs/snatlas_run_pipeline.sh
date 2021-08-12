#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=05:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N single_snuc_proc
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/snatlas_run_pipeline.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/snatlas_run_pipeline.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/run_singleprocessing.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/snrnaseq_data/" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/objects/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/initial_qc/" \
        --outs_structure TRUE;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/objects/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/initial_qc/";