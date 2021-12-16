#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -l mem=130gb
#PBS -S /bin/bash
#PBS -N spatial_tfs
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/estimate_spatial_tfs.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/estimate_spatial_tfs.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1.1_spatial_QC_reading/estimate_tfacts.R;