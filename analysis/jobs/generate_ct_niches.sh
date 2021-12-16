#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=05:00:00
#PBS -l mem=130gb
#PBS -S /bin/bash
#PBS -N annotate_obj
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ct_niches.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ct_niches.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/5_colocalization/find_niches_ct.R