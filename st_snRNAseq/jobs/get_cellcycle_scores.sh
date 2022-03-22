#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=30:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N cellcycle
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellcycle_scores.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellcycle_scores.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/estimate_cellcycle.R;
