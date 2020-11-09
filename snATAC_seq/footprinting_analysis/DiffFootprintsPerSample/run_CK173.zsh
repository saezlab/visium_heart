#!/usr/bin/env zsh
#
#### Job name
#SBATCH -J run_CK173
#SBATCH -e ./run_CK173.txt
#SBATCH -o ./run_CK173.txt
#SBATCH -t 30:00:00
#SBATCH --mem=180G -A rwth0429 -c 48


source ~/.zshrc
conda activate r-heart

split_bam=../SplitedBAMFiles/BAM
bam_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/BAM
bigwig_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/BigWig
peaks_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/Peaks
footprint_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/Footprints
motifmatching_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/MotifMatching
diff_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK173/Diff

#mkdir -p ${bigwig_loc}
#mkdir -p ${peaks_loc}
#mkdir -p ${footprint_loc}
#mkdir -p ${motifmatching_loc}
#mkdir -p ${bam_loc}
#mkdir -p ${diff_loc}
#
#samtools merge --threads 36 ${bam_loc}/Cardiomyocytes.bam ${split_bam}/CK173_0.bam ${split_bam}/CK173_1.bam ${split_bam}/CK173_2.bam
#samtools merge --threads 36 ${bam_loc}/Fibroblasts_1.bam ${split_bam}/CK173_3.bam ${split_bam}/CK173_5.bam
#samtools merge --threads 36 ${bam_loc}/Macrophages_1.bam ${split_bam}/CK173_4.bam
#samtools merge --threads 36 ${bam_loc}/Endothelial_1.bam ${split_bam}/CK173_6.bam
#samtools merge --threads 36 ${bam_loc}/Fibroblasts_2.bam ${split_bam}/CK173_7.bam
#samtools merge --threads 36 ${bam_loc}/Pericytes_2.bam ${split_bam}/CK173_8.bam
#samtools merge --threads 36 ${bam_loc}/Endothelial_2.bam ${split_bam}/CK173_9.bam
#samtools merge --threads 36 ${bam_loc}/Lymphatic_Endothelial.bam ${split_bam}/CK173_10.bam
#samtools merge --threads 36 ${bam_loc}/Neuronal.bam ${split_bam}/CK173_11.bam
#
#footprinting(){
#    	filename=$1
#    	samtools index $1
#    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g hs --nomodel -f BAMPE -q 0.05 --keep-dup all
#    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
#    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_summits.bed
#    	bamCoverage -bs 24 -b $filename -o $3/${filename%.bam}.bw -p 24 --normalizeUsing CPM
#	rgt-hint footprinting --organism=hg38 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
#	rgt-motifanalysis matching --organism=hg38 --output-location=$5 --input-files $4/${filename%.bam}.bed
#}
#
#cd ${bam_loc}
#for filename in *.bam
#do
#	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} 
#done

cd ${diff_loc}
rgt-hint differential --organism=hg38 --bc --nc 48 \
--mpbs-files=${motifmatching_loc}/Cardiomyocytes_mpbs.bed,\
${motifmatching_loc}/Fibroblasts_1_mpbs.bed,\
${motifmatching_loc}/Macrophages_1_mpbs.bed,\
${motifmatching_loc}/Endothelial_1_mpbs.bed,\
${motifmatching_loc}/Fibroblasts_2_mpbs.bed,\
${motifmatching_loc}/Pericytes_2_mpbs.bed,\
${motifmatching_loc}/Endothelial_2_mpbs.bed,\
${motifmatching_loc}/Lymphatic_Endothelial_mpbs.bed,\
${motifmatching_loc}/Neuronal_mpbs.bed \
--reads-files=${bam_loc}/Cardiomyocytes.bam,\
${bam_loc}/Fibroblasts_1.bam,\
${bam_loc}/Macrophages_1.bam,\
${bam_loc}/Endothelial_1.bam,\
${bam_loc}/Fibroblasts_2.bam,\
${bam_loc}/Pericytes_2.bam,\
${bam_loc}/Endothelial_2.bam,\
${bam_loc}/Lymphatic_Endothelial.bam,\
${bam_loc}/Neuronal.bam \
--conditions=Cardiomyocytes,Fibroblasts_1,Macrophages_1,Endothelial_1,Fibroblasts_2,Pericytes_2,Endothelial_2,Lymphatic_Endothelial,Neuronal \
--output-location=${diff_loc} \
--output-prefix=All





