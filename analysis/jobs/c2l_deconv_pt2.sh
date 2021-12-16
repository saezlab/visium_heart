#!/bin/bash
#SBATCH --job-name=nbsingularity
#SBATCH -t 2000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/jobs/c2l_deconv_pt2.out

module load system/singularity;

cd /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/;

for i in {AKK001_157785,AKK002_157779,AKK002_157781,AKK002_157782,AKK003_157775,AKK003_157777,AKK004_157772}
   do
   singularity exec --nv -B \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/ \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/cell2location-v0.05-alpha.sif \
   /bin/bash -c \
   "echo $i;
    python ./scripts/run_c2l.py $i";
   done
