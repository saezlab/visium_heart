#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N subsample_atlas
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/subsample_atlas.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/subsample_atlas.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/mini_tasks/make_atlas_subsamples.R;
