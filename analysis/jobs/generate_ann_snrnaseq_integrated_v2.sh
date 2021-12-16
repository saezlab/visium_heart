#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=05:00:00
#PBS -l mem=200gb
#PBS -S /bin/bash
#PBS -N annotate_obj
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ann_snrnaseq_integrated_v2.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ann_snrnaseq_integrated_v2.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/process_ann_object.R \
        --data_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann_umap.pdf";