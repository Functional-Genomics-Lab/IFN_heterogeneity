#!/bin/bash

# trim adapter sequences using cutadapt software
# adapter sequence is the Tn5 homology region in the 3' end of the read
# for pair-ended reads
SOURCE_DIR=~//interferon_project/ATAC-seq/fastq 

for i in ${SOURCE_DIR}/*_R1_001.fastq.gz
do
  SAMPLE=$(echo ${i} | sed "s/_R1_001\.fastq\.gz//")
  cutadapt -a GACAGAGAATATGTGTAGA -A TCTACACATATTCTCTGTC ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz -o ${SAMPLE}_trimmed_R1.fastq.gz -p ${SAMPLE}_trimmed_R2.fastq.gz -j 12  >${SAMPLE}_report_PE.txt 
done  
