#!/usr/bin/bash

source ~/.bashrc
conda activate r-4.0.3


sed '/var.2/d;/var.3/d' ../DiffFootprintsAllCellTypes/MotifMatching/${1}_mpbs.bed > ./${1}_mpbs.bed
python tf_to_gene.py ./${1}_mpbs.bed ../DiffFootprintsAllCellTypes//BAM/${1}.bam ./${1}.txt
rm ./${1}_mpbs.bed
