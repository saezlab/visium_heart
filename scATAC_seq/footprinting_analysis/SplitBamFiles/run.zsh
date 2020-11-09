#!/usr/local_rwth/bin/zsh

source ~/.zshrc
conda activate r-heart

python ./split_bam_by_cluster.py ./Clusters/$1.txt ../../run20_ATAC_visium/$1/outs/possorted_bam.bam ./BAM $1


