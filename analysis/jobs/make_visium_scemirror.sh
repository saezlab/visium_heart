#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -l mem=130gb
#PBS -S /bin/bash
#PBS -N scemirror
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_visium_scemirror.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_visium_scemirror.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/utils/sce_mirrors.R \
        --seurat_file "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds" \
        --sce_folder "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides_sce" \
        --assay "Spatial" \
        --reduction "umap_harmony";