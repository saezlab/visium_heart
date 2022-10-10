#!/bin/bash

### Job name
#SBATCH -J Neuronal
#SBATCH -e ./analysis.txt
#SBATCH -o ./analysis.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=100G -c 64 -x big-node001

source ~/.bashrc
conda activate r-4.1

#jupyter nbconvert --to html --execute ./01_subset.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./02_integrate_snATAC.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./03_integrate_snRNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./04_coembed_snATAC_snRNA.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./05_clustering.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./06_cleaning.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./07_annotation.ipynb --output-dir ../viz
jupyter nbconvert --to html --execute ./08_map_annotatiton_to_spatial.ipynb --output-dir ../viz

