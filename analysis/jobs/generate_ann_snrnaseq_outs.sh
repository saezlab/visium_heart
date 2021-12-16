#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=05:00:00
#PBS -l mem=90gb
#PBS -S /bin/bash
#PBS -N annotate_obj
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ann_snrnaseq_outs.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/generate_ann_snrnaseq_outs.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/annotate_object.R \
        --data_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples.rds" \
        --dictionary_path "/beegfs/work/hd_wh241/MI_revisions/markers/snrna_seq_cluster_annotations.txt" \
        --object_id "opt_clust_integrated" \
        --dictionary_id "opt_clust_integrated" \
        --new_variable "cell_type" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds";

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --used_assay "RNA" \
        --id_label "cell_type" \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/";

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/process_ann_object.R \
        --data_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann_umap.pdf";

$CONDA_PREFIX/bin/Rscript ./analysis/utils/h5ad_mirrors.R \
        --scell_data "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.h5ad";

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



