#!/bin/bash

#SBATCH --job-name=c2ltest
#SBATCH -t 1000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/c2l_nb_gpu.out

source ~/.bashrc;
conda activate cellpymc;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/;
$CONDA_PREFIX/bin/python nb_estimates_mjr.py;

