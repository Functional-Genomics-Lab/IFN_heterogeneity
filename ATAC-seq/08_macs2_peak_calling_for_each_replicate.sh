#!/bin/bash
# load modules
module load macs/2.1.2

# peak calling from aligned, filtered, duplicate removed seperate bed files (for each replicate, 14 in total)
SOURCE_DIR=~/project/Gozde_data/ATACseq/aligned_reads/Hg38
OUT_DIR=~/project/Gozde_data/ATACseq/peaks

for i in ${SOURCE_DIR}/*.bed
do

	SAMPLE=$(echo ${i} | sed "s/.*\///" |sed "s/\.bed//")

 
	macs2 callpeak -f BED --shift -75 --extsize 150 --cutoff-analysis --nomodel --call-summits --keep-dup all -p 0.01 -t $i -n $SAMPLE --tempdir ~/project/Gozde_data/ATACseq/ --outdir $OUT_DIR 2> $OUT_DIR/$SAMPLE"_macs2.log"
done  
