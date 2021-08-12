#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N snrnaseq_mrkrs
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_snrnaseq_mrkrplt.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/run_snrnaseq_mrkrplt.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --used_assay "RNA" \
        --id_label "opt_clust_integrated" \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples.rds" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/";