#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=18:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N snrnaseq_mrkr_estimation
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/estimate_snrseqclust_markers.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/estimate_snrseqclust_markers.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/3_functional_characterization/estimate_dea.R \
        --data_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples.rds" \
        --out_df "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_mrkrs.rds" \
        --group_class "opt_clust_integrated" \
        --test_assays "RNA" \
        --lfc "0.5" \
        --only_pos "yes";