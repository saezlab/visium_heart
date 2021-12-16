#!/bin/bash
#SBATCH --job-name=nbsingularity
#SBATCH -t 2000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/jobs/c2l_deconv_pt4.out

module load system/singularity;

cd /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/;

for i in {Visium_14_CK292,Visium_15_CK293,Visium_16_CK294,Visium_17_CK295,Visium_18_CK296,Visium_19_CK297,Visium_20_CK298}
   do
   singularity exec --nv -B \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/ \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/cell2location-v0.05-alpha.sif \
   /bin/bash -c \
   "echo $i;
    python ./scripts/run_c2l.py $i";
   done
