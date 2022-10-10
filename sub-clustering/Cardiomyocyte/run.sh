#!/bin/bash

### Job name
#SBATCH -J CM
#SBATCH -e ./analysis.txt
#SBATCH -o ./analysis.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=180G -c 64 -x big-node001

source ~/.bashrc
conda activate r-4.1

#jupyter nbconvert --to html --execute ./01_subset.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./02_integrate_snATAC.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./03_integrate_snRNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./04_coembed_snATAC_snRNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./05_clustering.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./06_map_cluster_to_spatial.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./07_cell_proportion.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./08_quality_check.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./09_cleaning.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./10_annotation.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./11_map_annotatiton_to_spatial.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./12_matching_ATAC_RNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./13_build_trajectory.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./14_subset_ATAC.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./15_tf_analysis.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./16_peak2gene.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./17_grn.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./18_function_analysis.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./19_grn.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./20_viz_TFs.ipynb --output-dir ../viz
jupyter nbconvert --to html --execute ./21_viz_TFs_in_spatial.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./22_viz_TFs.ipynb --output-dir ../viz
