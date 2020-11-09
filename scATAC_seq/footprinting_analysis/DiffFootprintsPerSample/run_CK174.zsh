#!/usr/bin/env zsh
#
#### Job name
#SBATCH -J run_CK174
#SBATCH -e ./run_CK174.txt
#SBATCH -o ./run_CK174.txt
#SBATCH -t 30:00:00
#SBATCH --mem=180G -c 48


source ~/.zshrc
conda activate r-heart

split_bam=../SplitedBAMFiles/BAM
bam_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/BAM
bigwig_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/BigWig
peaks_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/Peaks
footprint_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/Footprints
motifmatching_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/MotifMatching
diff_loc=/hpcwork/izkf/projects/SingleCellOpenChromatin/local/VisiumHeart/HINT/DiffFootprintsPerSample/CK174/Diff

mkdir -p ${bigwig_loc}
mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}

#samtools merge --threads 36 ${bam_loc}/Cardiomyocytes.bam ${split_bam}/CK174_0.bam 
#samtools merge --threads 36 ${bam_loc}/Fibroblasts_0.bam ${split_bam}/CK174_1.bam 
#samtools merge --threads 36 ${bam_loc}/Endothelial_1.bam ${split_bam}/CK174_2.bam 
#samtools merge --threads 36 ${bam_loc}/Pericytes.bam ${split_bam}/CK174_3.bam 

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
rgt-hint differential --organism=hg38 --bc --nc 100 \
--mpbs-files=${motifmatching_loc}/Cardiomyocytes_mpbs.bed,\
${motifmatching_loc}/Fibroblasts_0_mpbs.bed,\
${motifmatching_loc}/Endothelial_1_mpbs.bed \
--reads-files=${bam_loc}/Cardiomyocytes.bam,\
${bam_loc}/Fibroblasts_0.bam,\
${bam_loc}/Endothelial_1.bam \
--conditions=Cardiomyocytes,Fibroblasts_0,Endothelial_1 \
--output-location=${diff_loc} \
--output-prefix=All





