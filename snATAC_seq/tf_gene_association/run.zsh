#!/usr/bin/env zsh

source ~/.zshrc
conda activate r-heart

sample=$1
celltype=$2

if [ -f ../HINT/DiffFootprintsPerSample/${sample}/MotifMatching/${celltype}_mpbs.bed ];then
	# sed '/var.2/d;/var.3/d;/::/d' ../HINT/DiffFootprintsPerSample/${sample}/MotifMatching/${celltype}_mpbs.bed > ./${sample}/${celltype}_mpbs.bed
        sed '/var.2/d;/var.3/d' ../HINT/DiffFootprintsPerSample/${sample}/MotifMatching/${celltype}_mpbs.bed > ./${sample}/${celltype}_mpbs.bed
	python tf_to_gene.py ./${sample}/${celltype}_mpbs.bed ../HINT/DiffFootprintsPerSample/${sample}/BAM/${celltype}.bam ./${sample}/${celltype}.txt
	rm ./${sample}/${celltype}_mpbs.bed
fi
