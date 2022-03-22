#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=04:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N snrnapseudobulk
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_snrna_pseudobulk.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_snrna_pseudobulk.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --vars "cell_type" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds" \
        --def_assay "RNA";

$CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --vars "orig.ident,cell_type" \
        --collapsed "yes" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/psxsmpl_integrated_rnasamples_ann.rds" \
        --def_assay "RNA";
