#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=08:00:00
#PBS -l mem=150gb
#PBS -S /bin/bash
#PBS -N cell_states
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_Fib.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_Fib.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/utils/sce_mirrors.R \
        --seurat_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/Fib/Fib_states.rds" \
        --sce_folder "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/Fib/Fib_states_sce" \
        --assay "RNA" \
        --reduction "umap";

