#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=5:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N visium_mrkr_estimation
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/fetch_state_annotations.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/fetch_state_annotations.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/10_cellstates/fetch_annotations.R;