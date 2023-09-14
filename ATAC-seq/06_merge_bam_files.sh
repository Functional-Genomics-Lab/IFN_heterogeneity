# merge bam files which belong to the same condition (sample)
cd ~/project/Gozde_data/ATACseq/aligned_reads/Hg38
mkdir merged_filtered_bams 
#AI
samtools merge -@ 48 merged_filtered_bams/AI_merged_filtered.bam AI_210_S14_hg38_aligned_mtRemoved.bam AI_216_S11_hg38_aligned_mtRemoved.bam AI_220_S12_hg38_aligned_mtRemoved.bam AI_224_S13_hg38_aligned_mtRemoved.bam  
#AU
samtools merge -@ 48 merged_filtered_bams/AU_merged_filtered.bam AU_27_S10_hg38_aligned_mtRemoved.bam AU_215_S7_hg38_aligned_mtRemoved.bam AU_219_S8_hg38_aligned_mtRemoved.bam AU_223_S9_hg38_aligned_mtRemoved.bam
#DI
samtools merge -@ 48 merged_filtered_bams/DI_merged_filtered.bam DI_214_S4_hg38_aligned_mtRemoved.bam DI_217_S5_hg38_aligned_mtRemoved.bam DI_222_S6_hg38_aligned_mtRemoved.bam
#DU
samtools merge -@ 48 merged_filtered_bams/DU_merged_filtered.bam DU_23_S2_hg38_aligned_mtRemoved.bam DU_213_S1_hg38_aligned_mtRemoved.bam DU_221_S3_hg38_aligned_mtRemoved.bam 

