#!/bin/bash

#SBATCH --job-name=c2lrun
#SBATCH -t 1500:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/c2l_gpu.out

source ~/.bashrc;
conda activate cellpymc;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/;
$CONDA_PREFIX/bin/python run_c2l_ind.py;

