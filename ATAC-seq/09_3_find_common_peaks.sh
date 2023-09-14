#!/bin/bash
# load modules
module load bedtools
SOURCE_DIR=~/project/Gozde_data/ATACseq/peaks
# declare an array called array and define sample types
array=( AI AU DI DU )
SAMPLES=()
SAMPLES_BED=()

for j in ${array[@]}
do
	for i in $SOURCE_DIR/$j*_sorted.bed
	do
		#get sample ID for naming the new bed file
		SAMPLE_ID=$(echo $i | sed "s/.*peaks\\///" | sed "s/_sorted\.bed//")
		echo $SOURCE_DIR/merged_peaks_with_new_parameters/$j/*sorted.bed
		echo $i
		echo $SOURCE_DIR/${SAMPLE_ID}_intersected.bed
		# get only the peaks which are common in all AI, AU, etc. peaks
		bedtools intersect -wa -a $SOURCE_DIR/merged_peaks_with_new_parameters/$j/*sorted.bed -b ${i} > $SOURCE_DIR/${SAMPLE_ID}_intersected.bed 
		echo $SAMPLE_ID
	done
done 


