#!/bin/bash
# load modules 
module load bedtools

# sort bed files
sortBed -i ~/project/Gozde_data/ChIP_data/ENCFF710GQV.bed > ~/project/Gozde_data/ChIP_data/ENCFF710GQV_sorted.bed

sortBed -i ~/project/Gozde_data/ChIP_data/ENCFF361XMX.bed > ~/project/Gozde_data/ChIP_data/ENCFF361XMX_sorted.bed

sortBed -i ~/project/Gozde_data/ChIP_data/ENCFF695ETB.bed > ~/project/Gozde_data/ChIP_data/ENCFF695ETB_sorted.bed

sortBed -i ~/project/Gozde_data/ChIP_data/ENCFF725UFY.bed > ~/project/Gozde_data/ChIP_data/ENCFF725UFY_sorted.bed

sortBed -i ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_open.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_open_sorted.bed

sortBed -i ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_closed.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_closed_sorted.bed

sortBed -i ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_open.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_open_sorted.bed

sortBed -i ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_closed.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_closed_sorted.bed

# find the common peaks between histone  marks (H3K27Ac and H3K9me1)  which represent chromatin accessibility 
bedtools intersect -wa -wb -a ~/project/Gozde_data/ChIP_data/ENCFF710GQV_sorted.bed -b ~/project/Gozde_data/ChIP_data/ENCFF361XMX_sorted.bed > ~/project/Gozde_data/ChIP_data/regions_with_H3K27Ac_H3K4me1.bed

# do it for closed region histone marks (H3K27me3 and H3K9me3)
bedtools intersect -wa -wb -a ~/project/Gozde_data/ChIP_data/ENCFF695ETB_sorted.bed -b ~/project/Gozde_data/ChIP_data/ENCFF725UFY_sorted.bed > ~/project/Gozde_data/ChIP_data/regions_with_H3K27me3_H3K9me3.bed

# get peaks which is open in DAR however does not have open histone marks
bedtools subtract -A -a ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_open_sorted.bed -b ~/project/Gozde_data/ChIP_data/regions_with_H3K27Ac_H3K4me1.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_opened_withAPH.bed

bedtools subtract -A -a ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_open_sorted.bed -b ~/project/Gozde_data/ChIP_data/regions_with_H3K27Ac_H3K4me1.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_opened_withAPH.bed

bedtools subtract -A -a ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_closed_sorted.bed -b ~/project/Gozde_data/ChIP_data/regions_with_H3K27me3_H3K9me3.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_closed_withAPH.bed

bedtools subtract -A -a ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_closed_sorted.bed -b ~/project/Gozde_data/ChIP_data/regions_with_H3K27me3_H3K9me3.bed > ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_closed_withAPH.bed
