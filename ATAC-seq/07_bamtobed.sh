#!/bin/bash
# load modules
module load bedtools
# # convert bam files to bed
SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38

for i in ${SOURCE_DIR}/*_hg38_aligned_mtRemoved.bam
do
    SAMPLE_ID=$(echo ${i} |sed "s/.*Hg38\\///"| sed "s/_hg38_aligned_mtRemoved\.bam//")

    bedtools bamtobed -i ${i} > ${SOURCE_DIR}/${SAMPLE_ID}.bed
done
