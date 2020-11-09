#!/usr/bin/env zsh
#
#### Job name
#SBATCH -J run_CK170
#SBATCH -e ./run_CK170.txt
#SBATCH -o ./run_CK170.txt
#SBATCH -t 30:00:00
#SBATCH --mem=180G -A rwth0455 -c 48


source ~/.zshrc
conda activate r-heart

split_bam=../SplitedBAMFiles/BAM
bam_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/BAM
bigwig_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/BigWig
peaks_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/Peaks
footprint_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/Footprints
motifmatching_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/MotifMatching
diff_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK170/Diff

mkdir -p ${bigwig_loc}
mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}

#samtools merge --threads 36 ${bam_loc}/Endothelial_1.bam ${split_bam}/CK170_0.bam 
#samtools merge --threads 36 ${bam_loc}/Fibroblasts_3.bam ${split_bam}/CK170_1.bam ${split_bam}/CK170_5.bam
#samtools merge --threads 36 ${bam_loc}/Cardiomyocytes.bam ${split_bam}/CK170_2.bam 
#samtools merge --threads 36 ${bam_loc}/Pericytes.bam ${split_bam}/CK170_3.bam 
#samtools merge --threads 36 ${bam_loc}/Fibroblasts_5.bam ${split_bam}/CK170_4.bam 

footprinting(){
    	filename=$1
    	samtools index $1
    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g hs --nomodel -f BAMPE -q 0.05 --keep-dup all
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_summits.bed
    	bamCoverage -bs 24 -b $filename -o $3/${filename%.bam}.bw -p 24 --normalizeUsing CPM
	rgt-hint footprinting --organism=hg38 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
	rgt-motifanalysis matching --organism=hg38 --output-location=$5 --input-files $4/${filename%.bam}.bed
}

#cd ${bam_loc}
#for filename in *.bam
#do
#	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} 
#done

cd ${diff_loc}
rgt-hint differential --organism=hg38 --bc --nc 48 \
--mpbs-files=${motifmatching_loc}/Endothelial_1_mpbs.bed,\
${motifmatching_loc}/Fibroblasts_3_mpbs.bed,\
${motifmatching_loc}/Cardiomyocytes_mpbs.bed,\
${motifmatching_loc}/Fibroblasts_5_mpbs.bed \
--reads-files=${bam_loc}/Endothelial_1.bam,\
${bam_loc}/Fibroblasts_3.bam,\
${bam_loc}/Cardiomyocytes.bam,\
${bam_loc}/Fibroblasts_5.bam \
--conditions=Endothelial_1,Fibroblasts_3,Cardiomyocytes,Fibroblasts_5 \
--output-location=${diff_loc} \
--output-prefix=All





