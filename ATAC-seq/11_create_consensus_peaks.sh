#!/bin/bash
module load bedtools
# Merge all the intersected bed files to one consensus peak file:
SOURCE_DIR=~/project/Gozde_data/ATACseq/peaks

# Remove the column names
sed 1d $SOURCE_DIR/AI_common_peaks.bed > $SOURCE_DIR/AI_common_peaks2.bed
sed 1d $SOURCE_DIR/AU_common_peaks.bed > $SOURCE_DIR/AU_common_peaks2.bed
sed 1d $SOURCE_DIR/DI_common_peaks.bed > $SOURCE_DIR/DI_common_peaks2.bed
sed 1d $SOURCE_DIR/DU_common_peaks.bed > $SOURCE_DIR/DU_common_peaks2.bed

# concatenate bed files
cat  $SOURCE_DIR/AI_common_peaks2.bed $SOURCE_DIR/AU_common_peaks2.bed $SOURCE_DIR/DI_common_peaks2.bed $SOURCE_DIR/DU_common_peaks2.bed > $SOURCE_DIR/all_common_peaks.bed

# merge to get the consensus peaks
bedtools sort -i $SOURCE_DIR/all_common_peaks.bed > $SOURCE_DIR/all_common_peaks_sorted.bed
bedtools merge -i $SOURCE_DIR/all_common_peaks_sorted.bed > $SOURCE_DIR/consensus_peaks.bed

# Need to convert the bed file to saf format before creating the count matrix with featureCounts in R.
# create a new file with .saf extention which is basically another column with dots "."
# Don't forget to create the peak ID as the first column because you'll need the peak ID (Chr:Start-End) for featurecounts later!!!
awk 'OFS="\t" {print $1 ":" $2 "-" $3, $1, $2, $3, "."}' $SOURCE_DIR/consensus_peaks.bed > $SOURCE_DIR/consensus_peaks.saf
awk 'OFS="\t" {print $1 ":" $2 "-" $3, $1, $2, $3, "."}' $SOURCE_DIR/consensus_peaks_filt.bed > $SOURCE_DIR/consensus_peaks_filt.saf




