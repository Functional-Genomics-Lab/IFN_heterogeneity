#!/bin/bash
# load the appropreate modules
module load samtools                                                                      
# remove mt reads and low quality reads (MQ>30) from aligned files              
SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38
for i in ${SOURCE_DIR}/DU*_marked_duplicates.bam
                                            
do
  SAMPLE=$(echo ${i} | sed "s/_marked_duplicates\.bam//")
                                                                   
        # -v option in grep gets everything else except the given pattern       
        samtools view -h -@ 48 ${i} | grep -v 'chrM' | samtools view -@ 48 -b -h -F 4 -f 0x2 -q 30 -> ${SAMPLE}_aligned_mtRemoved.bam
done  
