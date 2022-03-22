#!/bin/bash
#SBATCH --job-name=nbsingularity
#SBATCH -t 2000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/jobs/c2l_deconv_pt3.out

module load system/singularity;

cd /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/;

for i in {Visium_7_CK285,Visium_8_CK286,Visium_9_CK287,Visium_10_CK288,Visium_11_CK289,Visium_12_CK290,Visium_13_CK291}
   do
   singularity exec --nv -B \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/ \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/cell2location-v0.05-alpha.sif \
   /bin/bash -c \
   "echo $i;
    python ./scripts/run_c2l.py $i";
   done
