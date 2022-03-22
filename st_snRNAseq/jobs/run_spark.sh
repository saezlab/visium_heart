#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=03:00:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N spatial_variability
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spark.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_spark.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1.1_spatial_QC_reading/run_spark.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_visium/objects/" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/results/spark/";