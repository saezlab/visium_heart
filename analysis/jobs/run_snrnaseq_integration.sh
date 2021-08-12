#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N snrnaseq_int
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_snrnaseq_integration.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_snrnaseq_integration.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/integrate_objects.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/objects/" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples.pdf" \
        --def_assay "RNA" \
        --default_resolution 1;