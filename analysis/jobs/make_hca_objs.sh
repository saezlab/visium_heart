#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N hcafiltering
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_hca_objs.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_hca_objs.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/mini_tasks/make_hca_objs.R;