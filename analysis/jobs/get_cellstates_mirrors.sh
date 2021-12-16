#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=05:00:00
#PBS -l mem=130gb
#PBS -S /bin/bash
#PBS -N cell_states
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_mirrors.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/get_cellstates_mirrors.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

for i in {cardiomyocyte,endothelial,fibroblast}
do
   $CONDA_PREFIX/bin/Rscript ./analysis/utils/h5ad_mirrors.R \
        --scell_data "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/${i}_states.rds" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/${i}_states.h5ad";

   $CONDA_PREFIX/bin/Rscript ./analysis/utils/sce_mirrors.R \
        --seurat_file "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/${i}_states.rds" \
        --sce_folder "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/${i}_states_sce" \
        --assay "RNA" \
        --reduction "umap";

   $CONDA_PREFIX/bin/Rscript ./analysis/6_atlas_comparison/get_pseudobulk_profiles.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/${i}_states.rds" \
        --vars "orig.ident,opt_state" \
        --collapsed "yes" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/results/ct_data/${i}/ps_${i}_states.rds" \
        --def_assay "RNA";

done
