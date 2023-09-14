##########################################################################
# plot Distribution of mapped reads
library(Rsamtools)
library(ggplot2)
library(magrittr)

bam_AI_210 <- "~/project/Gozde_data/ATACseq/aligned_reads/Hg38/AI_210_S14_hg38_aligned_mtRemoved.bam"

head(idxstatsBam(bam_AI_210))
#  seqnames seqlength  mapped unmapped
#1     chr1 248956422 6131337        0
#2     chr2 242193529 4796476        0
#3     chr3 198295559 3975231        0
#4     chr4 190214555 2646546        0
#5     chr5 181538259 3255076        0
#6     chr6 170805979 3698589        0

dim(idxstatsBam(bam_AI_210)) #[1] 195   4

# plot only the mapped read distribution to the important chromosomes
pdf("~/Gozde_data/ATACseq/files/AI_210_dist_mapped_reads.pdf", height= 15, width=28) 
idxstatsBam(bam_AI_210)[idxstatsBam(bam_AI_210)$seqnames %in% c(paste0("chr", 1:23), "chrX", "chrY"),] %>%
  ggplot(aes(seqnames,mapped,fill=seqnames))+
  geom_bar(stat="identity")+
  ylab("Number of mapped reads")+
  ggtitle("AI_210")+
  theme(plot.title = element_text(size = 40, face = "bold"),
        axis.text = element_text(size = 20), axis.title=element_text(size=25),
        axis.title.x = element_blank(),
        panel.border=element_rect(linetype=1,fill=NA),
        legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3))
dev.off()

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),fix="start",1)
TSSs

#####################################################################
###GenomicAlignments TSS enrichment
## AI_merged bam
library(GenomicAlignments)
library(ggplot2)
library(dplyr)

AI_merged = readGAlignmentPairs(
  "~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams/AI_merged_filtered.bam",
  param = ScanBamParam(
    mapqFilter = 1, 
    flag = scanBamFlag(
      isPaired = TRUE, 
      isProperPair = TRUE), 
    what = c("mapq", "isize")))


## AU_merged_bam
AU_merged = readGAlignmentPairs(
  "~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams/AU_merged_filtered.bam",
  param = ScanBamParam(
    mapqFilter = 1, 
    flag = scanBamFlag(
      isPaired = TRUE, 
      isProperPair = TRUE), 
    what = c("mapq", "isize")))


## DU_merged bam
DU_merged = readGAlignmentPairs(
  "~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams/DU_merged_filtered.bam",
  param = ScanBamParam(
    mapqFilter = 1, 
    flag = scanBamFlag(
      isPaired = TRUE, 
      isProperPair = TRUE), 
    what = c("mapq", "isize")))


## DI_merged bam
DI_merged = readGAlignmentPairs(
  "~/project/Gozde_data/ATACseq/aligned_reads/Hg38/merged_filtered_bams/DI_merged_filtered.bam",
  param = ScanBamParam(
    mapqFilter = 1, 
    flag = scanBamFlag(
      isPaired = TRUE, 
      isProperPair = TRUE), 
    what = c("mapq", "isize")))

################################################################
# generate TSS enrichment plot

pdf("~/project/Gozde_data/ATACseq/files/AI_merged_TSS_enrichment.pdf")
# put data in a dataframe to plot
df <- data.frame(x_axis_val=100*(-9:10-.5), values=tsse$values)
ggplot(df, aes(x_axis_val, values))+ geom_point(data = df, aes(y = values), colour = 'purple', size = 2) + geom_line(colour="purple") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  xlab("distance to TSS") +  ylab("aggregate TSS score")

dev.off() 


#########################################################################
### pie charts of functional elements
## AI_intersected peaks
# packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# data
AI_peaks <- read.table('~/project/Gozde_data/ATACseq/peaks/AI_intersected.bed')
colnames(AI_peaks) = c('chr','start','end','name','score','strand')
head(AI_peaks)
#chr  start    end              name score strand
#1 chr1 180795 180980 AI_210_S14_peak_1    41      .
#2 chr1 180795 180980 AI_210_S14_peak_1    41      .
#3 chr1 180800 180980 AI_210_S14_peak_1    41      .
#4 chr1 180795 180976 AI_210_S14_peak_1    41      .
#5 chr1 181354 181585 AI_210_S14_peak_2    54      .
#6 chr1 181354 181585 AI_210_S14_peak_2    54      .

(AI_peaks_gr = makeGRangesFromDataFrame(AI_peaks,keep.extra.columns=TRUE))
(MacsCalls_AI_filteredAnno <-  annotatePeak(AI_peaks_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene))

# plot
pdf("~/project/Gozde_data/ATACseq/files/Pie_chart_func_ele_AI_intersected_peaks.pdf")
plotAnnoPie(MacsCalls_AI_filteredAnno)
dev.off()
plotAnnoBar(MacsCalls_AI_filteredAnno) # bar plot form

## AU_intersected peaks
# packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# data
AU_peaks <- read.table('~/project/Gozde_data/ATACseq/peaks/AU_intersected.bed')
colnames(AU_peaks) = c('chr','start','end','name','score','strand')
head(AU_peaks)
#   chr  start    end             name score strand
#1 chr1  10076  10283 AU_215_S7_peak_1    87      .
#2 chr1  10078  10283 AU_215_S7_peak_1    87      .
#3 chr1  10084  10283 AU_215_S7_peak_1    87      .
#4 chr1 180756 181028 AU_215_S7_peak_2   120      .
#5 chr1 180756 181016 AU_215_S7_peak_2   120      .
#6 chr1 180799 181004 AU_215_S7_peak_2   120      .

(AU_peaks_gr = makeGRangesFromDataFrame(AU_peaks,keep.extra.columns=TRUE))
(MacsCalls_AU_filteredAnno <-  annotatePeak(AU_peaks_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene))

# plot
pdf("~/project/Gozde_data/ATACseq/files/Pie_chart_func_ele_AU_intersected_peaks.pdf")
plotAnnoPie(MacsCalls_AU_filteredAnno)
dev.off()
plotAnnoBar(MacsCalls_AU_filteredAnno) # bar plot form

## DI_intersected peaks
# packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# data
DI_peaks <- read.table('~/project/Gozde_data/ATACseq/peaks/DI_intersected.bed')
colnames(DI_peaks) = c('chr','start','end','name','score','strand')
head(DI_peaks)
#chr  start    end              name score strand
#1 chr1 10076 10404 DI_214_S4_peak_1   157      .
#2 chr1 10076 10404 DI_214_S4_peak_1   157      .
#3 chr1 10076 10404 DI_214_S4_peak_1   157      .
#4 chr1 10076 10404 DI_214_S4_peak_1   157      .
#5 chr1 17400 17582 DI_214_S4_peak_2   108      .
#6 chr1 17400 17582 DI_214_S4_peak_2   108      .

(DI_peaks_gr = makeGRangesFromDataFrame(DI_peaks,keep.extra.columns=TRUE))
(MacsCalls_DI_filteredAnno <-  annotatePeak(DI_peaks_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene))

# plot
pdf("~/project/Gozde_data/ATACseq/files/Pie_chart_func_ele_DI_intersected_peaks.pdf")
plotAnnoPie(MacsCalls_DI_filteredAnno)
dev.off()
plotAnnoBar(MacsCalls_DI_filteredAnno) # bar plot form

## DU_intersected peaks
# packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# data
DU_peaks <- read.table('~/project/Gozde_data/ATACseq/peaks/DU_intersected.bed')
colnames(DU_peaks) = c('chr','start','end','name','score','strand')
head(DU_peaks)
#chr  start    end              name score strand
#1 chr1 10072 10341 DU_213_S1_peak_1   172      .
#2 chr1 10072 10284 DU_213_S1_peak_1   172      .
#3 chr1 10072 10341 DU_213_S1_peak_1   172      .
#4 chr1 17399 17576 DU_213_S1_peak_2   106      .
#5 chr1 17408 17576 DU_213_S1_peak_2   106      .
#6 chr1 17408 17575 DU_213_S1_peak_2   106      .

(DU_peaks_gr = makeGRangesFromDataFrame(DU_peaks,keep.extra.columns=TRUE))
(MacsCalls_DU_filteredAnno <-  annotatePeak(DU_peaks_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene))

# plot
pdf("~/project/Gozde_data/ATACseq/files/Pie_chart_func_ele_DU_intersected_peaks.pdf")
plotAnnoPie(MacsCalls_DU_filteredAnno)
dev.off()
plotAnnoBar(MacsCalls_DU_filteredAnno) # bar plot form
###############################################################################
## PCA plot
# packages
library(DESeq2)
library(rtracklayer)
library(RColorBrewer)
library(ggplot2)

# data
load("~/project/Gozde_data/ATACseq/R_objs/myCounts.RData")

# define the groups
Group <- factor(c("AI","AI","AI","AI","AU","AU", "AU", "AU","DI", "DI", "DI", "DU","DU", "DU"))

# group replicates per condition
(metaData <- data.frame(Group, row.names = colnames(myCounts)))
#Group
#AI_210    AI
#AI_216    AI
#AI_220    AI
#AI_224    AI
#AU_215    AU
#AU_219    AU
#AU_223    AU
#AU_27     AU
#DI_214    DI
#DI_217    DI
#DI_222    DI
#DU_213    DU
#DU_221    DU
#DU_23     DU

# get the common peaks
consensus_peaks_filt <- read.table(file="~/project/Gozde_data/ATACseq/peaks/consensus_peaks_filt.bed")
# convert consensus_peaks_filt dataframe to GRanges object
colnames(consensus_peaks_filt) <- c("Chr", "Start", "End", "Width", "Strand")
consensus_peaks_filt_GR <- makeGRangesFromDataFrame(consensus_peaks_filt)

# create a DESeq2 object
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensus_peaks_filt)
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)

# save
save(atacDDS, file="~/project/Gozde_data/ATACseq/R_objs/atacDDS.RData")
# run it separately in terminal

# load
load(file="~/project/Gozde_data/ATACseq/R_objs/atacDDS.RData")

pca_atac <- plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))+
  scale_colour_brewer(palette = "Dark2", direction =-1, name = "Samples", labels = c("AI", "AU", "DI", "DU"),
                      guide = guide_legend(reverse=TRUE))+
  theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
        panel.border=element_rect(linetype=1,fill=NA),
        panel.grid.major = element_line(color = "grey",size = 0.1,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3))
ggsave(pca_atac, filename= "~/Gozde_data/ATACseq/files/PCA.pdf", width = 7, height = 12, device = cairo_pdf) 

############################################################################
#### Bar plots for GO term enrichment analysis
# packages
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggrastr)
library(dplyr)

### AUvsDU
## Open Regions
# load data
# don't get the first column as it's a duplicate of the second column
AUvsDU_open <- read.table(file="~/project/Gozde_data/ATACseq/files/APHMinusDMSO_open_GREAT_GO_BiolProcess.txt", sep = "\t")
dim(AUvsDU_open) #[1] 7316  13
c <- AUvsDU_open[1,] # get the colnames
AUvsDU_open <- AUvsDU_open[-1,] # remove the first col
dim(AUvsDU_open) #[1] 7315  13
head(AUvsDU_open)
colnames(AUvsDU_open) <- c
head(AUvsDU_open)

AUvsDU_open$minus_log10_padj <- -log10(as.numeric(AUvsDU_open$Hyper_Adjp_BH))

# sort based on -log10(padj) values
AUvsDU_open_ordered <- AUvsDU_open[order(AUvsDU_open$minus_log10_padj, decreasing=TRUE),]
# check 
head(AUvsDU_open_ordered)
AUvsDU_open_ordered$GeneRatio <- as.numeric(AUvsDU_open_ordered$Hyper_Foreground_Gene_Hits)/
  as.numeric(AUvsDU_open_ordered$Hyper_Background_Gene_Hits)

# round it to two decimals
AUvsDU_open_ordered$GeneRatio <- format(round(AUvsDU_open_ordered$GeneRatio,2), nsmall = 2)
AUvsDU_open_ordered$minus_log10_padj <- format(round(AUvsDU_open_ordered$minus_log10_padj,2), nsmall = 2)

pdf("~/Gozde_data/ATACseq/files/barplot_GO_BP_GM_AUvsDU_open.pdf", height= 10, width=18) 
AUvsDU_open_ordered[1:20,] %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  ggplot(aes(x=minus_log10_padj, y=name)) +
  geom_bar(aes(x=minus_log10_padj, y=name), stat="identity", fill= "red", color="black")+
  geom_text(aes(label = minus_log10_padj), color = "black",
            hjust = -0.1,
            size = 7,
            position = position_dodge(0.9))+
  scale_x_continuous(expand =c(0, 0),limits = c(0,40))+ # to fit all the values
  scale_y_discrete(expand = c(0, 0))+ # to remove space between y axis and origin
  ggtitle("Biological Process") + xlab("-log10(padj)")+ theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 40, face = "bold"),
        axis.text = element_text(size = 20), axis.title=element_text(size=25),
        axis.title.y = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
        legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3))

dev.off() 

### AUvsDU
## Closed Regions
# load data
# don't get the first column as it's a duplicate of the second column
AUvsDU_closed <- read.table(file="~/project/Gozde_data/ATACseq/files/APHMinusDMSO_closed_GREAT_GO_BiolProcess.txt", sep = "\t")
dim(AUvsDU_closed) #[1] 6954  13
c <- AUvsDU_closed[1,] # get the colnames
AUvsDU_closed <- AUvsDU_closed[-1,] # remove the first col
dim(AUvsDU_closed) #[1] 6953  13
head(AUvsDU_closed)
colnames(AUvsDU_closed) <- c
head(AUvsDU_closed)


AUvsDU_closed$minus_log10_padj <- -log10(as.numeric(AUvsDU_closed$Hyper_Adjp_BH))

# sort based on -log10(padj) values
AUvsDU_closed_ordered <- AUvsDU_closed[order(AUvsDU_closed$minus_log10_padj, decreasing=TRUE),]
# check 
head(AUvsDU_closed_ordered)
AUvsDU_closed_ordered$GeneRatio <- as.numeric(AUvsDU_closed_ordered$Hyper_Foreground_Gene_Hits)/
  as.numeric(AUvsDU_closed_ordered$Hyper_Background_Gene_Hits)

# round it to two decimals
AUvsDU_closed_ordered$GeneRatio <- format(round(AUvsDU_closed_ordered$GeneRatio,2), nsmall = 2)

pdf("~/Gozde_data/ATACseq/files/barplot_GO_BP_GM_AUvsDU_closed.pdf", height= 10, width=18) 
AUvsDU_closed_ordered[1:20,] %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  ggplot(aes(x=minus_log10_padj, y=name)) +
  geom_bar(aes(x=minus_log10_padj, y=name), stat="identity", fill= "blue", color="black")+
  geom_text(aes(label = minus_log10_padj), color = "black",
            hjust = -0.1,
            size = 7,
            position = position_dodge(0.9))+
  scale_x_continuous(expand =c(0, 0),limits = c(0,30))+ # to fit all the values
  scale_y_discrete(expand = c(0, 0))+ # to remove space between y axis and origin
  ggtitle("Biological Process") + xlab("-log10(padj)")+ theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 40, face = "bold"),
        axis.text = element_text(size = 20), axis.title=element_text(size=25),
        axis.title.y = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
        legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3))

dev.off() 
####################################################################
### Bar plots
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(reshape2)
## Plot how many DARs are in each comparison 
## Calculate the percentages of these numbers within the consensus peaks

##### Open DARs
### AUvsDU
# load data
load(file="~/project/Gozde_data/ATACseq/R_objs/APHMinusDMSO_U.RData")
load(file="~/project/Gozde_data/ATACseq/R_objs/APHMinusDMSO_I.RData")
load(file="~/project/Gozde_data/ATACseq/R_objs/AUvsDI.RData")
load(file="~/project/Gozde_data/ATACseq/R_objs/AIvsDU.RData")
consensus_peaks_filt <- read.table(file="~/project/Gozde_data/ATACseq/peaks/consensus_peaks_filt.bed")
nrow(consensus_peaks_filt) #[1] 203574

# check the number significantly accessible regions
length(APHMinusDMSO_U[(!is.na(APHMinusDMSO_U$padj) & 
                         APHMinusDMSO_U$padj < 0.05) & APHMinusDMSO_U$log2FoldChange > 0 ,]) # 21234
length(APHMinusDMSO_I[(!is.na(APHMinusDMSO_I$padj) & 
                         APHMinusDMSO_I$padj < 0.05)  & APHMinusDMSO_I$log2FoldChange > 0 ,]) #17993
length(AUvsDI[(!is.na(AUvsDI$padj) & AUvsDI$padj < 0.05) & AUvsDI$log2FoldChange > 0 ,]) #[1] 17586
length(AIvsDU[(!is.na(AIvsDU$padj) & AIvsDU$padj < 0.05) & AIvsDU$log2FoldChange > 0 ,]) #[1] 23286

# check the number significantly closed regions
length(APHMinusDMSO_U[(!is.na(APHMinusDMSO_U$padj) & 
                         APHMinusDMSO_U$padj < 0.05) & APHMinusDMSO_U$log2FoldChange < 0 ,]) # 30165
length(APHMinusDMSO_I[(!is.na(APHMinusDMSO_I$padj) & 
                         APHMinusDMSO_I$padj < 0.05) & APHMinusDMSO_I$log2FoldChange < 0 ,]) #20057
length(AUvsDI[(!is.na(AUvsDI$padj) & AUvsDI$padj < 0.05) & AUvsDI$log2FoldChange < 0,]) #[1] 24059
length(AIvsDU[(!is.na(AIvsDU$padj) & AIvsDU$padj < 0.05) & AIvsDU$log2FoldChange < 0,]) #[1] 26276

open_closed_dataf <- function(GR_name){
  library(ggplot2)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # annotate each DAR object
  open_sig <- GR_name[GR_name$padj < 0.05 & 
                               GR_name$log2FoldChange > 0,]
  anno_APHMinusDMSO_U_open <- annotatePeak(open_sig, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  a <- as.data.frame(anno_APHMinusDMSO_U_open)
  dim(a)
  table(a$annotation)
  # group the annotaiton categories to three groups: Promoter, intragenic and intergenic
  a$short_annotation <- ifelse(grepl("Promoter", a$annotation),  "Promoter", 
                                 ifelse(grepl(paste(c("UTR", "Exon", "Intron"), collapse='|'), a$annotation),  "Intragenic",
                                        ifelse(grepl(paste(c("Downstream", "Distal"), collapse='|'), a$annotation),  "Intergenic", FALSE)))
  table(a$short_annotation)
  #Intergenic Intragenic   Promoter 
  #7351       8285       5598 
  
  # create unique ids for each peak
  a$uniq_id <- paste0(a$seqnames, ":", a$start, "-", a$end)
  a$status <- rep("open", times=nrow(a))
  
  # Closed regions
  closed_sig <- GR_name[GR_name$padj < 0.05 & 
                               GR_name$log2FoldChange < 0,]
  anno_APHMinusDMSO_U_closed <- annotatePeak(closed_sig, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  b <- as.data.frame(anno_APHMinusDMSO_U_closed)
  dim(b) #[1] 30165    20
  table(b$annotation)
  # group the annotaiton categories to three groups: Promoter, intragenic and intergenic
  b$short_annotation <- ifelse(grepl("Promoter", b$annotation),  "Promoter", 
                               ifelse(grepl(paste(c("UTR", "Exon", "Intron"), collapse='|'), b$annotation),  "Intragenic",
                                      ifelse(grepl(paste(c("Downstream", "Distal"), collapse='|'), b$annotation),  "Intergenic", FALSE)))
  table(b$short_annotation)
  #Intergenic Intragenic   Promoter 
  #8338       16021       5806 
  
  # create unique ids for each peak
  b$uniq_id <- paste0(b$seqnames, ":", a$start, "-", a$end)
  b$status <- rep("closed", times=nrow(b))
  
  t <- rbind(cbind(a$uniq_id, as.numeric(a$log2FoldChange), as.numeric(a$padj), a$short_annotation, a$status), 
             cbind(b$uniq_id, as.numeric(b$log2FoldChange), as.numeric(b$padj), b$short_annotation, b$status))
  colnames(t) <- c("unique_id", "log2FoldChange", "padj", "short_annotation", "status")
  return(t) }

#APHMinusDMSO_I, AIvsDU, AUvsDI
x <- cbind(open_closed_dataf(APHMinusDMSO_U), rep("AUvsDU", nrow(open_closed_dataf(APHMinusDMSO_U))))
y <- cbind(open_closed_dataf(APHMinusDMSO_I), rep("AIvsDI", nrow(open_closed_dataf(APHMinusDMSO_I))))
z <- cbind(open_closed_dataf(AUvsDI), rep("AUvsDI", nrow(open_closed_dataf(AUvsDI))))
u <- cbind(open_closed_dataf(AIvsDU), rep("AIvsDU", nrow(open_closed_dataf(AIvsDU))))

df_peaks <- as.data.frame(rbind(x, y, z, u))
dim(df_peaks) #[1] 180656      6
head(df_peaks)
colnames(df_peaks) <- c("unique_id", "log2FoldChange", "padj", "short_annotation", "status", "comparison")


my_summ <- df_peaks %>% group_by(comparison) %>%
  summarise(
  n = n(),
  ratio = n/nrow(df_peaks),
  annotation = short_annotation
)
options(scipen = 999)
library(scales)
my_summ2 <- as.data.frame(table(my_summ$comparison))
colnames(my_summ2) <- c("comparison", "n")
# plot
pdf("~/project/Gozde_data/ATACseq/files/barPlot_count_of_sig_DARs.pdf")
my_summ2 %>% mutate(comparison = factor(comparison, levels = c("AIvsDU", "AUvsDU", "AUvsDI", "AIvsDI"))) %>% 
  ggplot(aes(x=comparison, y=n, fill= comparison)) +
  geom_bar(aes(fill= comparison), stat="identity")+
    ggtitle("") + xlab("Comparison")+ theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 40, face = "bold"),
        axis.text = element_text(size = 20), axis.title=element_text(size=25),
        axis.title.y = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
        panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3))

dev.off() 

