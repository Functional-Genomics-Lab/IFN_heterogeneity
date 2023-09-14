#!/bin/bash
# convert aligned sam to bam to save space

DEST_DIR=~/Interferon_project/ATAC-seq/10-13-22_data/analysis/files/aligned_reads
for i in *sam
do
	SAMPLE=$(echo ${i} | sed "s/\.sam//")

	samtools view -S -b ${i} > ${DEST_DIR}/${SAMPLE}.bam
	samtools sort -@ 8 ${DEST_DIR}/${SAMPLE}.bam > ${DEST_DIR}/${SAMPLE}_sorted.bam
done