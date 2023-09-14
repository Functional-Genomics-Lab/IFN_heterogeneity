#!/bin/bash

# align all trimmed fastqs using bwa-mem
SOURCE_DIR=~/interferon_project/ATAC-seq/fastq

# For correct indexing, gunzip the fna.gz file then run the following:
bwa index ~/genomes/bwa_approved_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

for i in ${SOURCE_DIR}/*_R1_001.fastq.gz
do
   SAMPLE=$(echo ${i} | sed "s/_R1_001\.fastq\.gz//")
   
   # Add the -t $SLURM_NTASKS_PER_NODE option for a fast bwa-mem mapping.
   bwa mem ~/genomes/bwa_approved_hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${SAMPLE}_trimmed_R1.fastq.gz ${SAMPLE}_trimmed_R2.fastq.gz  -t $SLURM_NTASKS_PER_NODE >  ${SAMPLE}_hg38_aligned_PE.sam

done
