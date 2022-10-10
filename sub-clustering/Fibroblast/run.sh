#!/bin/bash

### Job name
#SBATCH -J Fib 
#SBATCH -e ./analysis.txt
#SBATCH -o ./analysis.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=100G -c 64 

source ~/.bashrc
conda activate r-4.1

#jupyter nbconvert --to html --execute ./11_map_annotatiton_to_spatial.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./13_trajectory_human.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./14_matching_ATAC_RNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./15_subset_ATAC.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./16_add_trajectory.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./17_add_chromVAR.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./18_tf_analysis.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./19_peak2gene.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./20_grn.ipynb --output-dir ../viz
jupyter nbconvert --to html --execute ./21_viz_TFs.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./22_function_analysis.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./23_grn.ipynb --output-dir ../viz
