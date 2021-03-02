#!/bin/bash

#SBATCH --job-name=gputest
#SBATCH -t 100:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/jobs/test_gpu.out

source ~/.bashrc;
conda activate cellpymc;
cd /net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/;
$CONDA_PREFIX/bin/python test_gpu.py;

