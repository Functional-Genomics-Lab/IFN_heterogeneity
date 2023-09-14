#!/bin/bash                                                                    
# Load the modules
module load picard/2.10.3 
# use sorted bam files                                                         
SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38
          
for i in ${SOURCE_DIR}/*sorted.bam

do
        SAMPLE=$(echo ${i} | sed "s/\.sorted\.bam//")

        # remove duplicate reads from aligned and filtered files                
	java -jar $PICARD_DIR/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=${i}  O=${SAMPLE}_marked_duplicates.bam M=${SAMPLE}_marked_dup_metrics.xt
done
