#!/bin/bash
# load modules
module load samtools

SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams
OUT_DIR=~/project/Gozde_data/ATACseq/peaks
for i in ${SOURCE_DIR}/*_merged_filtered.bam
do
    SAMPLE=$(echo ${i} | sed "s/.*\///" |sed "s/_merged_filtered\.bam//")
    echo ${SAMPLE}
	echo ${OUT_DIR}/${SAMPLE}.SeqDepthNorm.bw
	
    # index the bam file
	samtools index ${i}
	
	# convert to bigwig file
	bamtoobamCoverage --bam ${i} -o ${OUT_DIR}/${SAMPLE}.SeqDepthNorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX --extendReads
done
