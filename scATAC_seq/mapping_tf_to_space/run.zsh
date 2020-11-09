#!/usr/bin/env zsh
#
#### Job name
#SBATCH -J tf_to_spatial
#SBATCH -e ./tf_to_spatial.txt
#SBATCH -o ./tf_to_spatial.txt
#SBATCH -t 10:00:00
#SBATCH --mem=180G -c 48

source ~/.zshrc
conda activate r_dorothea

Rscript tf_to_spatial.R
