#!/bin/bash
#
#### Job name
#SBATCH -J diff_footprinting
#SBATCH -e ./diff_footprinting.txt
#SBATCH -o ./diff_footprinting.txt
#SBATCH -t 120:00:00
#SBATCH --mem=960G --cpus-per-task=96

source ~/.bashrc
conda activate r-4.0.3

bam_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/BAM
bigwig_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/BigWig
peaks_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/Peaks
footprint_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/Footprints
motifmatching_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/MotifMatching
diff_loc=/data/scATA/SingleCellOpenChromatin/local/VisiumHeartRevision/snATAC/HINT/DiffFootprintsAllCellTypes/Diff

mkdir -p ${bigwig_loc}
mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}


for celltype in CM Endo Fib Lymphoid Myeloid Pericyte vSMCs;
do
samtools merge -f --threads 96 ${bam_loc}/${celltype}.bam ../SplitBamFilesByCellType/BAM/*_${celltype}.bam
done


footprinting(){
    	filename=$1
    	samtools index $1
    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g hs --nomodel -f BAMPE -q 0.01 --keep-dup all
    	sed -i '/random\|chrUn\|GL\|KI\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
    	bamCoverage -bs 24 -b $filename -o $3/${filename%.bam}.bw -p 24 --normalizeUsing CPM
     rgt-hint footprinting --organism=hg38 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
     rgt-motifanalysis matching --organism=hg38 --output-location=$5 --input-files $4/${filename%.bam}.bed
}

cd ${bam_loc}
for filename in *.bam;
do
	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} &
done

wait

cd ${diff_loc}
rgt-hint differential --organism=hg38 --bc --nc 96 \
--mpbs-files=${motifmatching_loc}/CM_mpbs.bed,\
${motifmatching_loc}/Endo_mpbs.bed,\
${motifmatching_loc}/Fib_mpbs.bed,\
${motifmatching_loc}/Lymphoid_mpbs.bed,\
${motifmatching_loc}/Myeloid_mpbs.bed,\
${motifmatching_loc}/Pericyte_mpbs.bed,\
${motifmatching_loc}/vSMCs_mpbs.bed \
--reads-files=${bam_loc}/CM.bam,\
${bam_loc}/Endo.bam,\
${bam_loc}/Fib.bam,\
${bam_loc}/Lymphoid.bam,\
${bam_loc}/Myeloid.bam,\
${bam_loc}/Pericyte.bam,\
${bam_loc}/vSMCs.bam \
--conditions=CM,Endo,Fib,Lymphoid,Myeloid,Pericyte,vSMCs \
--output-location=${diff_loc} \
--output-prefix=All
