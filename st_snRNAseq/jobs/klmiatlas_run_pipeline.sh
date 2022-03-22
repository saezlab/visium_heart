#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=08:00:00
#PBS -l mem=130gb
#PBS -S /bin/bash
#PBS -N single_snuc_proc
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/klmiatlas_run_pipeline.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/klmiatlas_run_pipeline.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/run_singleprocessing.R \
        --folder \
        --path "/beegfs/work/hd_wh241/Human_Acute_MI_Rafael/" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/processed_RNA/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/initial_qc/" \
        --outs_structure TRUE;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/processed_RNA/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/initial_qc/";

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/integrate_objects_kl.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/processed_RNA/" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/kl_miatlas_integrated_data.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/kl_miatlas_integrated_data_qcint.pdf" \
	--def_assay "RNA" \
	--default_resolution 1;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --id_label "opt_clust_integrated" \
        --path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/kl_miatlas_integrated_data.rds" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/kl_miatlas/";