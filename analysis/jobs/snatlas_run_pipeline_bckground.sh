#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -l mem=200gb
#PBS -S /bin/bash
#PBS -N single_snuc_proc
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/snatlas_run_pipeline_bckgrnd.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/snatlas_run_pipeline_bckgrnd.err
#PBS -q smp
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/run_singleprocessing_bckgrndcor.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/snrnaseq_data/" \
        --out_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/objects/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/initial_qc/" \
        --outs_structure TRUE;

$CONDA_PREFIX/bin/Rscript ./analysis/1_snuc_QC_doublets_majorannotation/plot_knownmarkers.R \
        --folder \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/objects/" \
        --out_fig_path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/initial_qc/";

$CONDA_PREFIX/bin/Rscript ./analysis/2_snuc_integration_harmony/integrate_objects.R \
        --path "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/objects/" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/integration/integrated_rnasamples.rds" \
        --out_fig_file "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq_bckgrnd/integration/integrated_rnasamples.pdf" \
        --def_assay "RNA" \
        --batch_file "/beegfs/work/hd_wh241/MI_revisions/markers/snrna_batch_ann.csv" \
        --default_resolution 1;