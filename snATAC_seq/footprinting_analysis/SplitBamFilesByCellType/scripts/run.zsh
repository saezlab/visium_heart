#!/bin/sh

source ~/.bashrc
conda activate r-4.0.3

cluster_file=./$1.txt
bam_file=../../../Alignment/$1/outs/possorted_bam.bam
output_location=../BAM

python ~/scopen/scripts/split_bam_by_cluster.py $cluster_file $bam_file $output_location $1
