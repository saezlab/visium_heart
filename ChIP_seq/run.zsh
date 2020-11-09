#!/usr/bin/env zsh

source ~/.zshrc
conda activate r-heart

bowtie2_index=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/Bowtie2_indexes/hg38/hg38

sra=$1

# download SRA
prefetch $sra

# Dump each read into separate file. Files will receive suffix corresponding to read number.
fastq-dump -I --split-files $sra
rm -rf $sra
	
# Adapter sequences were trimmed from FASTQs
trim_galore --suppress_warn --cores 8 --paired -o ./ ${sra}_1.fastq ${sra}_2.fastq
rm ${sra}_1.fastq
rm ${sra}_2.fastq

# Map the reads to reference genome
bowtie2 --very-sensitive --no-discordant -x ${bowtie2_index} -1 ${sra}_1_val_1.fq -2 ${sra}_2_val_2.fq -S ${sra}.map.sam -X 2000 -p 50
rm ${sra}_1_val_1.fq
rm ${sra}_2_val_2.fq

# Filter out reads mapped to chrY, mitochondria, and unassembled "random" contigs, 
sed -i '/chrY/d;/chrM/d;/random/d;/chrUn/d' ${sra}.map.sam

# Convert sam file to bam file, sort the result and generate the index file
samtools view -Sb ${sra}.map.sam > ${sra}.map.bam
samtools sort ${sra}.map.bam -o ${sra}.sort.bam
samtools index ${sra}.sort.bam
rm ${sra}.map.sam
rm ${sra}.map.bam

# Remove duplicates 
java -jar /hpcwork/izkf/jar/picard.jar MarkDuplicates INPUT=${sra}.sort.bam \
OUTPUT=${sra}.rmdup.bam METRICS_FILE=${sra}_matrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
rm ${sra}.sort.bam
rm ${sra}.sort.bam.bai

# and bad map quality
samtools view -bq 30 ${sra}.rmdup.bam > ${sra}.filter.bam
samtools index ${sra}.filter.bam
rm ${sra}.rmdup.bam

# Require reads to be properly paired
samtools view -f2 ${sra}.filter.bam -b > ${sra}.bam
samtools index ${sra}.bam
samtools flagstat ${sra}.bam > ${sra}_qc.txt

rm ${sra}.filter.bam
rm ${sra}.filter.bam.bai
