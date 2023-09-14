#!/bin/bash
# load modules
module load bedtools

# merge all the peak files prior to DiffBind
SOURCE_DIR=~/project/Gozde_data/ATACseq/peaks

# declare an array called array and define sample types
array=( AI AU DI DU )
SAMPLES=()
SAMPLES_BED=()

for j in ${array[@]}
do
	# for i in $SOURCE_DIR/$j*.narrowPeak #code for each replicate peak file
	for i in $SOURCE_DIR/merged_peaks_with_new_parameters/$j/*.narrowPeak
	do
	#get sample ID for naming the new bed file
	#SAMPLE_ID=$(echo $i | sed "s/.*peaks\///" | sed "s/_peaks\.narrowPeak//") #code for each replicate peak file
	
#	echo $SAMPLE_ID
	# convert all the narroPeak files to bed files. Narrowpeak is in bed 6+4 format
	# so take only the first 6 cols
#	cut -f 1-6 ${i} > $SOURCE_DIR/$SAMPLE_ID.bed #code for each replicate peak file
	cut -f 1-6 ${i} > $SOURCE_DIR/merged_peaks_with_new_parameters/$j/${j}.bed
	
	# sort the newly formed bed files
#	bedtools sort -i $SOURCE_DIR/$SAMPLE_ID.bed > $SOURCE_DIR/${SAMPLE_ID}_sorted.bed #code for each replicate peak file
	bedtools sort -i $SOURCE_DIR/merged_peaks_with_new_parameters/$j/${j}.bed > $SOURCE_DIR/merged_peaks_with_new_parameters/$j/${j}_sorted.bed
#	SAMPLES+=( "$SOURCE_DIR/${SAMPLE_ID}_sorted.bed" ) #code for each replicate peak file
	echo ${i}
	done
echo ${j}
echo ${SAMPLES[@]}
echo ${SAMPLES[0]}

SAMPLES=()
done 

