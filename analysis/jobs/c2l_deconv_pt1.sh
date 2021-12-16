#!/bin/bash
#SBATCH --job-name=nbsingularity
#SBATCH -t 2000:00
#SBATCH --mail-user=roramirezf@uni-heidelberg.de
#SBATCH --mail-type=END
#SBATCH --output /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/jobs/c2l_deconv_pt1.out

module load system/singularity;

cd /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/;

for i in {AKK006_157771,Visium_1_CK279,Visium_2_CK280,Visium_3_CK281,Visium_4_CK282,Visium_5_CK283,Visium_6_CK284}
   do
   singularity exec --nv -B \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/ \
   /net/data.isilon/ag-saez/bq_rramirez/MI_deconvolution/cell2location-v0.05-alpha.sif \
   /bin/bash -c \
   "echo $i;
    python ./scripts/run_c2l.py $i";
   done
