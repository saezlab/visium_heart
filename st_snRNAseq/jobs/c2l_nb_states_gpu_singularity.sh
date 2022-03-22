#!/bin/bash

#SBATCH --job-name=nbsingularity
#SBATCH -t 1000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/jobs/c2l_nb_states_gpu_singularity.out

module load system/singularity;

cd /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/;

singularity exec --nv -B \
/net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/ \
/net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/cell2location-v0.05-alpha.sif \
/bin/bash -c \
"python ./scripts/nb_estimates_states_singularity.py";
