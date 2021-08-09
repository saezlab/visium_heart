#!/bin/bash

#SBATCH --job-name=getstates
#SBATCH -t 400:00
#SBATCH --mem=64GB
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/get_states.out

source ~/.bashrc;
conda activate scell_analysis;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions;
for i in {cardiomyocytes,endothelial_cells,fibroblasts,macrophages,pericytes}
do
   $CONDA_PREFIX/bin/Rscript get_cell_states.R \
        --data_path /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_annotated.rds \
        --out_path /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.rds \
        --out_fig_file /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.pdf \
        --cell_class $i \
        --class_label "cell_type"
        --start_res 0.2;

   $CONDA_PREFIX/bin/Rscript add_funcomics.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.rds" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.rds" \
        --group_class "opt_state" \
        --gene_set_collection "not";

    mkdir ./results/ct_data/${i};

    $CONDA_PREFIX/bin/Rscript estimate_dea.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.rds" \
        --out_df "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}/${i}_dea.rds" \
        --test_assays "RNA,progeny,dorothea" \
        --group_class "opt_state" \
        --lfc "0.20" \
        --only_pos "yes";

    $CONDA_PREFIX/bin/Rscript plot_funcomics.R \
        --data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}_states.rds" \
        --dea_data_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}/${i}_dea.rds" \
        --nfeats 10 \
        --ngenes_ORA 10 \
        --pvalue_ORA 0.15 \
        --pvalue 0.001 \
        --lfc 0.5 \
        --test_assays "RNA,progeny,dorothea" \
        --group_class "opt_state" \
        --out_path "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/ct_data/${i}/${i}";
done
