#!/bin/bash
# guide: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
module load homer/4.9

# AUvsDU significantly open regions (p_adj <0.05)
findMotifsGenome.pl ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_U_open.bed hg38 ~/project/Gozde_data/ATACseq/files/ -size 200 -bg ~/project/Gozde_data/ATACseq/peaks/consensus_peaks_filt.bed


# AIvsDI significantly open regions (p_adj <0.05)
findMotifsGenome.pl ~/project/Gozde_data/ATACseq/peaks/APHMinusDMSO_I_open.bed hg38 ~/project/Gozde_data/ATACseq/files/ -size 200 -bg ~/project/Gozde_data/ATACseq/peaks/consensus_peaks_filt.bed

# Possible enhancer regions near IFNs
findMotifsGenome.pl ~/project/Gozde_data/ATACseq/files/regions_of_interest/possible_enhancers_near_IFNs.bed hg38 ~/project/Gozde_data/ATACseq/files/regions_of_interest/motif_enrichment/ -size 200 -bg ~/project/Gozde_data/ATACseq/peaks/consensus_peaks_filt.bed
