#!/bin/bash

### Job name
#SBATCH -J Pericyte
#SBATCH -e ./analysis.txt
#SBATCH -o ./analysis.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 120:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=900G -c 64 -x big-node001

source ~/.bashrc
conda activate r-4.1

jupyter nbconvert --to html --execute ./01_integrate_snRNA.ipynb --output-dir ../viz
jupyter nbconvert --to html --execute ./02_clustering.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./06_annotation.ipynb --output-dir ../viz
#jupyter nbconvert --to html --execute ./07_map_annotatiton_to_spatial.ipynb --output-dir ../viz

