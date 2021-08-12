#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=01:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N spatial_pseudobulk
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_umap.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_umap.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/mini_tasks/get_umap_from_integration.R;