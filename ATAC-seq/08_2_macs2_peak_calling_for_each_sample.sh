#!/bin/bash
# load modules
module load macs/2.1.2

# peak calling from aligned, filtered, duplicate removed bam files
SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams
OUTDIR=~/project/Gozde_data/ATACseq/peaks/merged_peaks_with_new_parameters

for i in ${SOURCE_DIR}/*.bed
do

	SAMPLE=$(echo ${i} |sed "s/.*merged_filtered_bams\\///" |sed "s/\.bed//")
    macs2 callpeak --cutoff-analysis -f BED --shift -75 --extsize 150 --nomodel --call-summits --keep-dup all -p 0.01 -t ${i} -n ${SAMPLE} --outdir ${OUTDIR}/${SAMPLE}/ 2> ${OUTDIR}/${SAMPLE}/${SAMPLE}_macs2.log

done
