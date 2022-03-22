#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=12gb
#SBATCH --job-name=human_mi
#SBATCH --output=singularity_job-%j.out
#SBATCH --export=NONE

## This script will process all Fibroblast-Myeloid RNAscope images from Christoph Kuppe for spot counting and nuclear segmentation

## Define parameters
module load system/singularity
module load devel/java_jdk/1.8.0

singluarity_cache="/home/hd/hd_hd/hd_gr294/singularity_cache"
bfconvert_bin="/home/hd/hd_hd/hd_gr294/bin/bftools/bfconvert"

project_dir="/gpfs/bwfor/work/ws/hd_gr294/"
orig_image=$project_dir"/original_images"
single_channels=$project_dir"/single_channels"
mesmer_out=$project_dir"/mesmer"
nuclei_tables=$project_dir"/nuclei_tables"

declare -a channelarray=("0" "1" "2" "3")

cd $orig_image

## Run pipeline
for FILE in *.nd2
do
    file_base=(${FILE//./ })
    file_name=${file_base[0]}
    outfile_name=$file_name"nuclei_params.csv"
    if [ ! -f $nuclei_tables/$outfile_name ]
    then
        echo $FILE
        for channel in ${channelarray[@]}
        do
            echo $channel
            outfile_name_tif=$file_name".ch_"$channel".tif"
	    if [ ! -f $single_channels/$outfile_name_tif ]
	    then
	    	$bfconvert_bin -channel $channel $orig_image/$FILE $single_channels/$outfile_name_tif
	    fi
	done

	## Run mesmer nuclear detection
	singularity exec -B $single_channels:/input -B $mesmer_out:/output $singluarity_cache/vanvalenlab-deepcell-applications-0.3.0.img python /usr/src/app/run_app.py mesmer --nuclear-image /input/$file_name".ch_0".tif --nuclear-channel 0 --output-directory /output --output-name $file_name".mesmer_nuclear_mask.tif" --compartment nuclear --image-mpp 0.1 --squeeze

	## Get centroid positions of nuclei from mesmser masks
	python_script_dir="/gpfs/bwfor/work/ws/hd_gr294/scripts"
	singularity exec -B $mesmer_out:/input -B $nuclei_tables:/output -B $python_script_dir:/scripts $singluarity_cache/jupyter.scipy_notebook_hub.2_1_1.sif python /scripts/get_nuclei_params.py /input/$file_name".mesmer_nuclear_mask.tif" /output/$file_name"nuclei_params.csv"
	else
		echo $FILE" is already processed!"
	fi
done
