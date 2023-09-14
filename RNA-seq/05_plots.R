#### Create plots for the paper ####
# packages
library(edgeR)
library(DESeq2)
#BiocManager::install("reshape2")
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer) 
library(Cairo) 
#BiocManager::install("ggrepel")
library(ggrepel)
#BiocManager::install("ggrastr")
library(ggrastr)
library(pheatmap)
library("stringr")
#####################################################
# set the colors for samples
hex_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
#####################################################
# PCA plot
# GM12878
# get the data
load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_counts_keep.RData")
# create group names for GM samples
Group_GM <- factor(c("GM_DU","GM_DI", "GM_AU", "GM_AI", "GM_DU", "GM_DI", "GM_AU", "GM_AI", "GM_DU", "GM_DI", "GM_AU", "GM_AI"))
colnames(GM_counts_keep) <- Group_GM
load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_RNA_DDS.RData")

#normalize the reads
RNA_Rlog_GM <- rlog(GM_RNA_DDS)

# plot
(pca_gm <- plotPCA(RNA_Rlog_GM, intgroup = "Group_GM", ntop = nrow(RNA_Rlog_GM))+
  scale_colour_brewer(palette = "Dark2", direction =-1, name = "Samples", labels = c("AI", "AU", "DI", "DU"),
                       guide = guide_legend(reverse=TRUE))+
  theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
        panel.border=element_rect(linetype=1,fill=NA),
        panel.grid.major = element_line(color = "grey",size = 0.1,linetype = 2),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA), 
        axis.line = element_line(color = "black",size = 0.3)))
ggsave(pca_gm, filename= "~/Gozde_data/RNA-seq/files/PCA_GM.pdf", device = cairo_pdf) # 

  
# JURKAT
# get the data
load("~/project/Gozde_data/RNA-seq/R_objs/JU_counts_keep.RData")
# create group names for JU samples
Group_JU <- factor(c("JU_DU","JU_DI", "JU_AU", "JU_AI", "JU_DU", "JU_DI", "JU_AU", "JU_DI", "JU_AI", "JU_AU", "JU_AI", "JU_DU"))
colnames(JU_counts_keep) <- Group_JU
load(file="~/project/Gozde_data/RNA-seq/R_objs/JU_RNA_DDS.RData")

#normalize the reads
RNA_Rlog_JU <- rlog(JU_RNA_DDS)

# plot
(pca_ju <- plotPCA(RNA_Rlog_JU, intgroup = "Group_JU", ntop = nrow(RNA_Rlog_JU))+
    scale_colour_brewer(palette = "Dark2", direction =  -1,name = "Samples", labels = c("AI", "AU", "DI", "DU"),
                        guide = guide_legend(reverse=TRUE))+
    theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
          panel.border=element_rect(linetype=1,fill=NA),
          panel.grid.major = element_line(color = "grey",size = 0.1,linetype = 2),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA), 
          axis.line = element_line(color = "black",size = 0.3)))
ggsave(pca_ju, filename= "~/Gozde_data/RNA-seq/files/PCA_JU.pdf", device = cairo_pdf) # 

#################################################################################################
# GM12878
# box plot
# IFNB1
pdf("~/Gozde_data/RNA-seq/files/Boxplot_IFNB1_GM.pdf") # 
(p_IFNB1 <- ggplot(GM_df_TC[GM_df_TC$id=="IFNB1",], 
                   mapping =aes(x=sample, y=value)) + geom_boxplot(aes(x=sample, y=value, fill=sample),
                                                                   color="black")+
    scale_colour_brewer(palette = "Dark2")+
   geom_point()+
    ggtitle("IFNB1 in GM") + ylab("log2(CPM+1)")+
    scale_x_discrete(labels=c('DU', 'DI', 'AU', 'AI'))+
    theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
          axis.title.x = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
          legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA), 
          axis.line = element_line(color = "black",size = 0.3))+
    scale_y_continuous(limits = c(0,6)))
dev.off() 

# IFNL1
pdf("~/Gozde_data/RNA-seq/files/Boxplot_IFNL1_GM.pdf") # 
(p_IFNL1 <- ggplot(GM_df_TC[GM_df_TC$id=="IFNL1",], mapping =aes(x=sample, y=value)) + geom_boxplot(aes(x=sample, y=value, fill=sample), color="black")+
    scale_colour_brewer(palette = "Dark2")+
    geom_point()+
    ggtitle("IFNL1 in GM") + ylab("log2(CPM+1)")+
    scale_x_discrete(labels=c('DU', 'DI', 'AU', 'AI'))+
    theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
          axis.title.x = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
          legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA), 
          axis.line = element_line(color = "black",size = 0.3))+
    scale_y_continuous(limits = c(0,7.5)))
 dev.off() 

###############################################################
# JURKAT
JU_pseudo_TC <- log2(cpm(JU_counts_keep) + 1)

JU_df_TC <- melt(JU_pseudo_TC, id = rownames(JU_counts_keep))
names(JU_df_TC)[1:2] <- c ("id", "sample")
JU_df_TC$method <- rep("TC", nrow(JU_df_TC))

# box plot
# IFNB1
pdf("~/Gozde_data/RNA-seq/files/Boxplot_IFNB1_JU.pdf") 
(p_IFNB1 <- ggplot(JU_df_TC[JU_df_TC$id=="IFNB1",], mapping =aes(x=sample, y=value)) + geom_boxplot(aes(x=sample, y=value, fill=sample), color="black")+
   scale_colour_brewer(palette = "Dark2")+  
   ggtitle("IFNB1 in JU") + ylab("log2(CPM+1)")+
    geom_point()+
    scale_x_discrete(labels=c('DU', 'DI', 'AU', 'AI'))+
   theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
         axis.title.x = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
         legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
         panel.background = element_rect(fill = "transparent",colour = NA),
         plot.background = element_rect(fill = "transparent",colour = NA), 
         axis.line = element_line(color = "black",size = 0.3))+
    scale_y_continuous(limits = c(0,4.5)))
dev.off() 

# IFNL1
pdf("~/Gozde_data/RNA-seq/files/Boxplot_IFNL1_JU.pdf") 
(p_IFNL1 <- ggplot(JU_df_TC[JU_df_TC$id=="IFNL1",], mapping =aes(x=sample, y=value)) + geom_boxplot(aes(x=sample, y=value, fill=sample), color="black")+
   scale_colour_brewer(palette = "Dark2")+
    geom_point()+
   ggtitle("IFNL1 in JU") + ylab("log2(CPM+1)")+
    scale_x_discrete(labels=c('DU', 'DI', 'AU', 'AI'))+
   theme(axis.text = element_text(size = 20), axis.title=element_text(size=25),
         axis.title.x = element_blank(), panel.border=element_rect(linetype=1,fill=NA),
         legend.position = "none", panel.grid.major = element_line(color = "grey",size = 0.3,linetype = 2),
         panel.background = element_rect(fill = "transparent",colour = NA),
         plot.background = element_rect(fill = "transparent",colour = NA), 
         axis.line = element_line(color = "black",size = 0.3))+
    scale_y_continuous(limits = c(0,4.5)))
dev.off()



######################################################
rm(list = ls())
library(dplyr)
library(ggpubr)
library(Seurat)
library(reshape2)
library(destiny)
library(ggrastr)
library(DESeq2)
library(edgeR)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

# Load dataset
load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_RNA_DDS.RData")

# Normalize the reads
GM_RNA_Rlog <- rlog(GM_RNA_DDS)

# Extract results - AI-DU
GM_AI_DU = results(GM_RNA_DDS, c("Group_GM", "GM_AI", "GM_DU"), format = "DataFrame")
GM_AI_DU = cbind(Gene = rownames(GM_AI_DU), GM_AI_DU)

# Save results
rio::export(GM_AI_DU, file="~/project/Gozde_data/RNA-seq/files/GM_APH_IMinusDMSO_U.xlsx")


# Plot up-down genes in all comparisons
GM_AI_DU = results(GM_RNA_DDS, c("Group_GM", "GM_AI", "GM_DU"), format = "DataFrame")
GM_AU_DU = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
GM_DI_DU = results(GM_RNA_DDS, c("Group_GM", "GM_DI", "GM_DU"), format = "DataFrame")

GM_AI_DU = GM_AI_DU[!(is.na(GM_AI_DU$padj)),]
GM_AU_DU = GM_AU_DU[!(is.na(GM_AU_DU$padj)),]
GM_DI_DU = GM_DI_DU[!(is.na(GM_DI_DU$padj)),]

GM_AI_DU = cbind(Gene = rownames(GM_AI_DU), GM_AI_DU)
GM_AU_DU = cbind(Gene = rownames(GM_AU_DU), GM_AU_DU)
GM_DI_DU = cbind(Gene = rownames(GM_DI_DU), GM_DI_DU)

# Plot up-down genes in all
au_up = GM_AU_DU[GM_AU_DU$log2FoldChange > 0.6 & GM_AU_DU$padj < 0.05, 'Gene']
ai_up = GM_AI_DU[GM_AI_DU$log2FoldChange > 0.6 & GM_AI_DU$padj < 0.05,'Gene']
di_up = GM_DI_DU[GM_DI_DU$log2FoldChange > 0.6 & GM_DI_DU$padj < 0.05,'Gene']

au_down = GM_AU_DU[GM_AU_DU$log2FoldChange < -0.6 & GM_AU_DU$padj < 0.05,'Gene']
ai_down = GM_AI_DU[GM_AI_DU$log2FoldChange < -0.6 & GM_AI_DU$padj < 0.05,'Gene']
di_down = GM_DI_DU[GM_DI_DU$log2FoldChange < -0.6 & GM_DI_DU$padj < 0.05,'Gene']

toplot = data.frame(vars = c(rep('UP', 3), rep('DOWN', 3)),
		vals = c(length(au_up), length(ai_up), length(di_up), length(au_down), length(ai_down), length(di_down)),
			comp = rep(c('AU_DU', 'AI_DU', 'DI_DU'), 2))
toplot$comp = factor(toplot$comp, levels = c('AU_DU', 'DI_DU', 'AI_DU'))

pdf('GM_ALLCOMPs_UPDOWN.pdf', width = 8)
ggbarplot(toplot, x = 'comp', y = 'vals', color = 'vars', fill = 'vars', palette = c('red', 'blue'), position = position_dodge(0.9))+
	ylab('Number of genes') + xlab('') +
	theme(text=element_text(size=25, face = 'bold')) 
dev.off()


# Overlaps - AI-DU and AU-DU
library(VennDiagram)

pdf('GM_AIDU_AUDU_OVERLAP.pdf')
draw.pairwise.venn(area1 = length(c(ai_up, ai_down)),
			area2 = length(c(au_up, au_down)),
			cross.area = sum(sum(ai_up %in% au_up), sum(ai_down %in% au_down)),
			fill = c('orange', 'darkgreen'),
			col = c('orange', 'darkgreen'),
			category = c('AIDU', 'AUDU'), fontface = 'bold', cex = 2,
			cat.cex = 2, cat.fontface = 'bold', margin = 0.13)
dev.off()

pdf('GM_AIDU_DIDU_OVERLAP.pdf')
draw.pairwise.venn(area1 = length(c(ai_up, ai_down)),
			area2 = length(c(di_up, di_down)),
			cross.area = sum(sum(ai_up %in% di_up), sum(ai_down %in% di_down)),
			fill = c('orange', 'darkgreen'),
			col = c('orange', 'darkgreen'),
			category = c('AIDU', 'DIDU'), fontface = 'bold', cex = 2,
			cat.cex = 2, cat.fontface = 'bold', margin = 0.13)
dev.off()




# Plot up-down genes in AU-DU
GM_AU_DU <- results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")

vals = c(sum(GM_AU_DU$log2FoldChange > 0.6 & GM_AU_DU$padj < 0.05),
	sum(GM_AU_DU$log2FoldChange < -0.6 & GM_AU_DU$padj < 0.05))

toplot = data.frame(vars = c('UP', 'DOWN'), vals = vals)

pdf('GM_AU_DU_UPDOWN.pdf')
ggbarplot(toplot, x = 'vars', y = 'vals', color = 'vars', fill = 'vars', palette = c('blue', 'red')) +
	ylab('Number of genes') + xlab('') +
	theme(text=element_text(size=25, face = 'bold')) +
	NoLegend()
dev.off()


# Gene Ontology Enrichment AU-DU UP
au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames
bcg = rownames(GM_AI_DU)

gores_audu = GOenrich(au_du_up, bcg, qCut = 0.05)
gores_audu$ObsOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$ObsAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$BcgOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$BcgAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$OddsRatio = (gores_audu$ObsOv / gores_audu$ObsAll) / (gores_audu$BcgOv / gores_audu$BcgAll)

rio::export(gores_audu, "~/project/Gozde_data/RNA-seq/files/GO_Enrichment_AU_HIGHER_DU.xlsx")
gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/GO_Enrichment_AU_HIGHER_DU.xlsx")

gores_audu[gores_audu$GO == 'MF', 'Description'] %>% head
gores_audu[gores_audu$GO == 'BP', 'Description'] %>% head
gores_audu[gores_audu$GO == 'CC', 'Description'] %>% head

gonames = c("cytokine activity", "inflammatory response", "cytokine-mediated signaling pathway", "chemotaxis")
toplot = gores_audu[gores_audu$Description %in% gonames,]
toplot$log10FDR = -log10(toplot$qvalue)
toplot$Description = factor(toplot$Description, levels = rev(c(toplot$Description)))

pdf('GO_Enrichment_AU_UP_DU.pdf', width = 12)
ggscatter(toplot, x = 'log10FDR', y = 'Description', color = 'blue', size = 'OddsRatio', xlim = c(0, 25)) +
	ylab('') + xlab(expression('-log'[10]*'(FDR)')) +
	geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
	scale_size_continuous(range = c(6,10)) + 
	theme(text=element_text(size=25, face = 'bold'), legend.position = 'top')
dev.off()



# Gene Ontology Enrichment AU-DU DOWN
au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du_down = au_du[au_du$log2FoldChange < -0.6 & au_du$padj < 0.05,] %>% rownames
bcg = rownames(GM_AU_DU)

gores_audu = GOenrich(au_du_down, bcg, qCut = 0.05)
gores_audu$ObsOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$ObsAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$BcgOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$BcgAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$OddsRatio = (gores_audu$ObsOv / gores_audu$ObsAll) / (gores_audu$BcgOv / gores_audu$BcgAll)

rio::export(gores_audu, "~/project/Gozde_data/RNA-seq/files/GO_Enrichment_AU_LOWER_DU.xlsx")
gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/GO_Enrichment_AU_LOWER_DU.xlsx")

gores_audu[gores_audu$GO == 'MF', 'Description'] %>% head
gores_audu[gores_audu$GO == 'BP', 'Description'] %>% head
gores_audu[gores_audu$GO == 'CC', 'Description'] %>% head

gonames = c("ribonucleoprotein complex binding", "ribosome biogenesis", "rRNA processing", "ribonucleoprotein complex biogenesis")
toplot = gores_audu[gores_audu$Description %in% gonames,]
toplot$log10FDR = -log10(toplot$qvalue)
toplot = toplot[order(toplot$log10FDR, decreasing = T),]
toplot$Description = factor(toplot$Description, levels = rev(c(toplot$Description)))

pdf('GO_Enrichment_AU_DOWN_DU.pdf', width = 12)
ggscatter(toplot, x = 'log10FDR', y = 'Description', color = 'red', size = 'OddsRatio', xlim = c(0, 40)) +
	ylab('') + xlab(expression('-log'[10]*'(FDR)')) +
	geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
	scale_size_continuous(range = c(6,10)) + 
	theme(text=element_text(size=25, face = 'bold'), legend.position = 'top')
dev.off()

## HEATMAPS ##
load("~/project/Gozde_data/RNA-seq/R_objs/GM_counts_keep.RData")

GM_counts_norm = log2(cpm(GM_counts_keep) + 1)
GM_counts_norm_noai <- GM_counts_norm[, !(grepl(c("_AI"), colnames(GM_counts_norm)))]

GM_counts_scaled_noai <- scale(t(GM_counts_norm_noai))
GM_counts_scaled_noai <- t(GM_counts_scaled_noai)

# Significantly upregulated genes
au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames

au_di = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DI"), format = "DataFrame")
au_di_up = au_di[au_di$log2FoldChange > 0.6 & au_di$padj < 0.05,] %>% rownames

# Chemokine ligand genes #
chemogns = read.csv('~/project/Gozde_data/RNA-seq/files/Chemokine_Ligands_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

#gns_toplot <- intersect(chemogns, rownames(GM_counts_scaled_noai))
gns_toplot <- intersect(intersect(chemogns, au_du_up), au_di_up)

mat_plotgns <- GM_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3), paste0('GM_DI_', 1:3),
										paste0('GM_AU_', 1:3)))

pdf("heatMap_GM_DUandAU_CCL_CXCLs_SIGNIFICANT.pdf", width = 9)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()

# IL genes #

# Interleukin gene symbols in human
ilgns = read.csv('~/project/Gozde_data/RNA-seq/files/Interleukins_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

#gns_toplot <- intersect(ilgns, rownames(GM_counts_scaled_noai))
gns_toplot <- intersect(intersect(ilgns, au_du_up), au_di_up)

mat_plotgns <- GM_counts_scaled_noai[gns_toplot,]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3), paste0('GM_DI_', 1:3),
										paste0('GM_AU_', 1:3)))

pdf("heatMap_GM_DUandAU_Interleukins_SIGN.pdf", width = 9)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()



# TNF genes #

# Interleukin gene symbols in human
tnfgns = read.csv('~/project/Gozde_data/RNA-seq/files/TumorNecrosisFactors_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

#gns_toplot <- intersect(tnfgns, rownames(GM_counts_scaled_noai))
gns_toplot <- intersect(intersect(tnfgns, au_du_up), au_di_up)
mat_plotgns <- GM_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3), paste0('GM_DI_', 1:3),
										paste0('GM_AU_', 1:3)))

pdf("heatMap_GM_DUandAU_TNFs_SIGN.pdf", width = 9, height = 5)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()


# Cell Cycle Genes #
GM_counts_norm = log2(cpm(GM_counts_keep) + 1)
GM_counts_norm_audu <- GM_counts_norm[, !(grepl(c("_AI|_DI"), colnames(GM_counts_norm)))]

GM_counts_scaled_audu <- scale(t(GM_counts_norm_audu))
GM_counts_scaled_audu <- t(GM_counts_scaled_audu)

s_gns = intersect(cc.genes$s.genes, rownames(GM_counts_scaled_audu))
g2m_gns = intersect(cc.genes$g2m.genes, rownames(GM_counts_scaled_audu))

au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
s_gns_up = au_du[rownames(au_du) %in% s_gns & au_du$padj < 0.05 & au_du$log2FoldChange > 0, ] %>% rownames
s_gns_down = au_du[rownames(au_du) %in% s_gns & au_du$padj < 0.05 & au_du$log2FoldChange < 0, ] %>% rownames

g2m_gns_up = au_du[rownames(au_du) %in% g2m_gns & au_du$padj < 0.05 & au_du$log2FoldChange > 0, ] %>% rownames
g2m_gns_down = au_du[rownames(au_du) %in% g2m_gns & au_du$padj < 0.05 & au_du$log2FoldChange < 0, ] %>% rownames

# S
mat_plotgns <- GM_counts_scaled_audu[c(s_gns_up, s_gns_down),]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3), paste0('GM_AU_', 1:3)))

pdf("heatMap_S_Phase_Sign.pdf", width = 8, height = 10)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()

mat_plotgns <- GM_counts_norm_audu[c('CDC45', 'CDC6', 'PCNA', 'E2F8'),]
melted_int <- data.frame(melt(mat_plotgns))
melted_int$Group = gsub('_[0-9]', '', melted_int$Var2) %>% gsub('GM_', '', .)
melted_int$Group = factor(melted_int$Group, levels = c('DU', 'AU'))

pdf("heatMap_S_Phase_Sign_BoxPlot.pdf", width = 6, height = 8)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="Normalized Gene Expression")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  facet_wrap(~Var1, scales = 'free') + NoLegend() +
  rotate_x_text(45)
dev.off()


# G2M
mat_plotgns <- GM_counts_scaled_audu[c(g2m_gns_up, g2m_gns_down),]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3), paste0('GM_AU_', 1:3)))

pdf("heatMap_G2M_Phase_Sign.pdf", width = 8, height = 13)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()





# Motif Enrichment Genes - Only AU and DU #
GM_counts_norm = log2(cpm(GM_counts_keep) + 1)
GM_counts_norm_noai <- GM_counts_norm[, !(grepl(c("_AI|_DI"), colnames(GM_counts_norm)))]
GM_counts_scaled_noai <- scale(t(GM_counts_norm_noai))
GM_counts_scaled_noai <- t(GM_counts_scaled_noai)

nfkb = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/nfkb_subunit_genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
fos = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/FOS_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
jun = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/JUN_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
bzip = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/bZIP_domain_genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
irfs = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/IRF_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

gns_toplot = c(nfkb, fos, jun, bzip, irfs)
gns_toplot <- intersect(gns_toplot, rownames(GM_counts_scaled_noai))

mat_plotgns <- GM_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('GM_DU_', 1:3),
										paste0('GM_AU_', 1:3)))

pdf("heatMap_GM_DUandAU_MOTIFGENES_ALL.pdf", width = 9, height = 15)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()


# Boxplots of motif enrichment genes
au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames

GM_counts_norm = log2(cpm(GM_counts_keep) + 1)

gns_toplot <- intersect(gns_toplot, au_du_up)
# plot only the genes significantly upregulated in AUvsDU in both GM and JU
gns_toplot <- c("FOS", "JUN", "ATF3", "CEBPB", "JDP2", "MAFF", "IRF1", "IRF7", "IRF9")
mat_plotgns <- GM_counts_norm[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
melted_int$Group = gsub('_[0-9]', '', melted_int$Var2) %>% gsub('GM_', '', .)

melted_int$Group = factor(melted_int$Group, levels = c('DU', 'DI', 'AU', 'AI'))
# change the level names for the paper
levels(melted_int$Group) <- c('CU', 'CI', 'AU', 'AI')

pdf("BoxPlot_MotifGenes.pdf", width = 9, height = 12)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold'),legend.position = "none") +
  facet_wrap(~Var1, scales = 'free') + 
  rotate_x_text(45)
dev.off()

# plot the common motifs between GM and JU
pdf("~/project/Gozde_data/RNA-seq/files/BoxPlot_MotifGenes_common_GM.pdf", width = 9, height = 12)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold'),legend.position = "none") +
  facet_wrap(~Var1, scales = 'free') + 
  rotate_x_text(45)
dev.off()

#############################################################################
## GM vs JU

# Load datasets
load(file="~/project/Gozde_data/RNA-seq/R_objs/JU_counts_keep.RData")
load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_counts_keep.RData")

# Create matrix that only has the DU samples
GM_counts_du = GM_counts_keep[, grepl('DU', colnames(GM_counts_keep))]
JU_counts_du = JU_counts_keep[, grepl('DU', colnames(JU_counts_keep))]

# Merge the matrices
hmm = merge(GM_counts_du, JU_counts_du, by="row.names", all=TRUE)
hmm[is.na(hmm)] = 0
rownames(hmm) = hmm[, 'Row.names']
hmm$Row.names = NULL
mergedMat = hmm

# create group names for GM samples
Groups = factor(c("GM_DU","GM_DU", "GM_DU", "JU_DU", "JU_DU", "JU_DU"))

# create DESeq2 object with filtered count matrix
metaData <- data.frame(Groups, row.names = colnames(mergedMat))
Merged_RNA_DDS <- DESeqDataSetFromMatrix(mergedMat, metaData, ~Groups)
Merged_RNA_DDS <- DESeq(Merged_RNA_DDS)

saveRDS(Merged_RNA_DDS, file="~/project/Gozde_data/RNA-seq/files/DU_GM_vs_JU_RNA_DDS.RDS")

# Normalize the reads
Merged_RNA_Rlog <- rlog(Merged_RNA_DDS)

# Extract results - AI-DU
DU_GM_JU = results(Merged_RNA_DDS, c("Groups", "GM_DU", "JU_DU"), format = "DataFrame")
DU_GM_JU = DU_GM_JU[!(is.na(DU_GM_JU$log2FoldChange)),]
#DU_GM_JU = cbind(Gene = rownames(DU_GM_JU), DU_GM_JU)

gm_ju_up = DU_GM_JU[DU_GM_JU$log2FoldChange > 1 & DU_GM_JU$padj < 0.05,] %>% rownames
gm_ju_down = DU_GM_JU[DU_GM_JU$log2FoldChange < -1 & DU_GM_JU$padj < 0.05,] %>% rownames

rio::export(as.data.frame(gm_ju_up), file="gm_ju_up.xlsx")
rio::export(as.data.frame(gm_ju_down), file="gm_ju_down.xlsx")

# Save results
rio::export(DU_GM_JU, file="~/project/Gozde_data/RNA-seq/files/DU_GM_vs_JU.xlsx")


#############################################################################
## GM and JU markers in AU vs DU


load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_RNA_DDS.RData")
load(file="~/project/Gozde_data/RNA-seq/R_objs/JU_RNA_DDS.RData")

gm_au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
gm_au_du = gm_au_du[!(is.na(gm_au_du$log2FoldChange)),]
gm_au_du_up = gm_au_du[gm_au_du$log2FoldChange > 0.6 & gm_au_du$padj < 0.05,] %>% rownames

ju_au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
ju_au_du = ju_au_du[!(is.na(ju_au_du$padj)),]
ju_au_du_up = ju_au_du[ju_au_du$log2FoldChange > 0.6 & ju_au_du$padj < 0.05,] %>% rownames

#length(intersect(gm_ju_up, ju_au_du_up)) / length(union(gm_ju_up, ju_au_du_up))
#length(intersect(gm_ju_down, ju_au_du_up)) / length(union(gm_ju_down, ju_au_du_up))

#length(intersect(gm_ju_up, gm_au_du_up)) / length(union(gm_ju_up, gm_au_du_up))
#length(intersect(gm_ju_down, gm_au_du_up)) / length(union(gm_ju_down, gm_au_du_up))


library(GeneOverlap)
gm_ju_list = list(gm_ju_up, gm_ju_down)
au_du_list = list(gm_au_du_up, ju_au_du_up)
names(gm_ju_list) = c('GM_UP_JU', 'GM_DOWN_JU')
names(au_du_list) = c('GM_AU_UP_DU', 'JU_AU_UP_DU')

#bcg = length(Reduce(union, list(rownames(DU_GM_JU), rownames(gm_au_du), rownames(ju_au_du))))
bcg = length(Reduce(intersect, list(rownames(DU_GM_JU), rownames(gm_au_du), rownames(ju_au_du))))

resgom = newGOM(gm_ju_list, au_du_list, genome.size=bcg)
oddsr = getMatrix(resgom, name="odds.ratio")
pvals = getMatrix(resgom, name="pval")

oddsrM = melt(oddsr)
pvalsM = melt(pvals)
pvalsM$FDR = p.adjust(as.numeric(pvalsM$value))

toplot = cbind(oddsrM, FDR = pvalsM$FDR)
toplot$value = ifelse(toplot$value > 1000, 1000, toplot$value)
toplot$is_top = ifelse(toplot$value > 3 & toplot$FDR > 3, 'Top', 'Others')
#toplot$Var1 = rownames(toplot)
toplot$FDR = formatC(toplot$FDR, format = "e", digits = 2)
toplot$log10FDR = -log10(p.adjust(as.numeric(pvalsM$value)))
toplot$Var2 = factor(toplot$Var2, levels = rev(c('OPC-OL_OPC', 'OPC-OL_Middle', 'OPC-OL_Oli', 'OPC-AST_OPC', 'OPC-AST_Middle', 'OPC-AST_Ast', 'MIC-AST_Ast', 'MIC-AST_Middle', 'MIC-AST_Mic')))

pdf("GMJU_vs_AUDU_Enrichment.pdf", width = 8)
ggplot(toplot, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white") +
 rotate_x_text(90) +
 theme_classic() +
 scale_fill_gradient2(midpoint = 1, high = 'red', low = 'white') +
 #geom_text(aes(label = round(oddsr, digits = 2)), fontface = 'bold') +
 geom_text(aes(label = FDR), fontface = 'bold') +
 ylab('') + xlab('') +
 theme(text=element_text(size=20, face = 'bold')) +
 #labs(fill = expression('-log'[10]*'(FDR)')) +
 labs(fill = 'Odds\nRatio') +
 coord_flip() + rotate_x_text(45)
dev.off()


#############################################################################
## JU Plot UP-DOWN Genes


# Plot up-down genes in all comparisons
JU_AI_DU = results(JU_RNA_DDS, c("Group_JU", "JU_AI", "JU_DU"), format = "DataFrame")
JU_AU_DU = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
JU_DI_DU = results(JU_RNA_DDS, c("Group_JU", "JU_DI", "JU_DU"), format = "DataFrame")

JU_AI_DU = JU_AI_DU[!(is.na(JU_AI_DU$padj)),]
JU_AU_DU = JU_AU_DU[!(is.na(JU_AU_DU$padj)),]
JU_DI_DU = JU_DI_DU[!(is.na(JU_DI_DU$padj)),]

JU_AI_DU = cbind(Gene = rownames(JU_AI_DU), JU_AI_DU)
JU_AU_DU = cbind(Gene = rownames(JU_AU_DU), JU_AU_DU)
JU_DI_DU = cbind(Gene = rownames(JU_DI_DU), JU_DI_DU)

# Plot up-down genes in all
au_up = JU_AU_DU[JU_AU_DU$log2FoldChange > 0.6 & JU_AU_DU$padj < 0.05, 'Gene']
ai_up = JU_AI_DU[JU_AI_DU$log2FoldChange > 0.6 & JU_AI_DU$padj < 0.05,'Gene']
di_up = JU_DI_DU[JU_DI_DU$log2FoldChange > 0.6 & JU_DI_DU$padj < 0.05,'Gene']

au_down = JU_AU_DU[JU_AU_DU$log2FoldChange < -0.6 & JU_AU_DU$padj < 0.05,'Gene']
ai_down = JU_AI_DU[JU_AI_DU$log2FoldChange < -0.6 & JU_AI_DU$padj < 0.05,'Gene']
di_down = JU_DI_DU[JU_DI_DU$log2FoldChange < -0.6 & JU_DI_DU$padj < 0.05,'Gene']

toplot = data.frame(vars = c(rep('UP', 3), rep('DOWN', 3)),
		vals = c(length(au_up), length(ai_up), length(di_up), length(au_down), length(ai_down), length(di_down)),
			comp = rep(c('AU_DU', 'AI_DU', 'DI_DU'), 2))
toplot$comp = factor(toplot$comp, levels = c('AU_DU', 'DI_DU', 'AI_DU'))

pdf('JU_ALLCOMPs_UPDOWN.pdf', width = 8)
ggbarplot(toplot, x = 'comp', y = 'vals', color = 'vars', fill = 'vars', palette = c('red', 'blue'), position = position_dodge(0.9))+
	ylab('Number of genes') + xlab('') +
	theme(text=element_text(size=25, face = 'bold')) 
dev.off()


# Overlaps - AI-DU and AU-DU
library(VennDiagram)

pdf('JU_AIDU_AUDU_OVERLAP.pdf')
draw.pairwise.venn(area1 = length(c(ai_up, ai_down)),
			area2 = length(c(au_up, au_down)),
			cross.area = sum(sum(ai_up %in% au_up), sum(ai_down %in% au_down)),
			fill = c('orange', 'darkgreen'),
			col = c('orange', 'darkgreen'),
			category = c('AIDU', 'AUDU'), fontface = 'bold', cex = 2,
			cat.cex = 2, cat.fontface = 'bold', margin = 0.13)
dev.off()

pdf('JU_AIDU_DIDU_OVERLAP.pdf')
draw.pairwise.venn(area1 = length(c(ai_up, ai_down)),
			area2 = length(c(di_up, di_down)),
			cross.area = sum(sum(ai_up %in% di_up), sum(ai_down %in% di_down)),
			fill = c('orange', 'darkgreen'),
			col = c('orange', 'darkgreen'),
			category = c('AIDU', 'DIDU'), fontface = 'bold', cex = 2,
			cat.cex = 2, cat.fontface = 'bold', margin = 0.13)
dev.off()

# Overlap with GM
load(file="~/project/Gozde_data/RNA-seq/R_objs/GM_RNA_DDS.RData")
GM_AU_DU = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
GM_AU_DU = GM_AU_DU[!(is.na(GM_AU_DU$padj)),]
GM_AU_DU = cbind(Gene = rownames(GM_AU_DU), GM_AU_DU)
au_up_gm = GM_AU_DU[GM_AU_DU$log2FoldChange > 0.6 & GM_AU_DU$padj < 0.05, 'Gene']
au_down_gm = GM_AU_DU[GM_AU_DU$log2FoldChange < -0.6 & GM_AU_DU$padj < 0.05,'Gene']

pdf('JU_GM_AUDU_OVERLAP.pdf')
draw.pairwise.venn(area1 = length(c(au_up, au_down)),
			area2 = length(c(au_up_gm, au_down_gm)),
			cross.area = sum(sum(au_up %in% au_up_gm), sum(au_down %in% au_down_gm)),
			fill = c('orange', 'darkgreen'),
			col = c('orange', 'darkgreen'),
			category = c('AUDU_JU', 'AUDU_GM'), fontface = 'bold', cex = 2,
			cat.cex = 2, cat.fontface = 'bold', margin = 0.13)
dev.off()

library(GeneOverlap)
resgom = newGOM(list(c(au_up, au_down)), list(c(au_up_gm, au_down_gm)), genome.size = nrow(GM_AU_DU))
pvalmat = getMatrix(resgom, name="pval")



# Plot up-down genes in AU-DU
JU_AU_DU <- results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
JU_AU_DU = JU_AU_DU[!(is.na(JU_AU_DU$padj)),]

vals = c(sum(JU_AU_DU$log2FoldChange > 0.6 & JU_AU_DU$padj < 0.05),
	sum(JU_AU_DU$log2FoldChange < -0.6 & JU_AU_DU$padj < 0.05))

toplot = data.frame(vars = c('UP', 'DOWN'), vals = vals)

pdf('JU_AU_DU_UPDOWN.pdf')
ggbarplot(toplot, x = 'vars', y = 'vals', color = 'vars', fill = 'vars', palette = c('blue', 'red')) +
	ylab('Number of genes') + xlab('') +
	theme(text=element_text(size=25, face = 'bold')) +
	NoLegend()
dev.off()


#############################################################################
## JU Gene Ontology Enrichments


# Gene Ontology Enrichment AU-DU UP
au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames
bcg = rownames(au_du)

gores_audu = GOenrich(au_du_up, bcg, qCut = 0.05)
gores_audu$ObsOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$ObsAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$BcgOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$BcgAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$OddsRatio = (gores_audu$ObsOv / gores_audu$ObsAll) / (gores_audu$BcgOv / gores_audu$BcgAll)

rio::export(gores_audu, "~/project/Gozde_data/RNA-seq/files/JU_GO_Enrichment_AU_HIGHER_DU.xlsx")
ju_gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/JU_GO_Enrichment_AU_HIGHER_DU.xlsx")
gm_gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/GM_GO_Enrichment_AU_HIGHER_DU.xlsx")

ju_gores_audu$CellType = 'Jurkat'
gm_gores_audu$CellType = 'GM'

# Examine top enrichments
gores_audu[gores_audu$GO == 'MF', 'Description'] %>% head
gores_audu[gores_audu$GO == 'BP', 'Description'] %>% head
gores_audu[gores_audu$GO == 'CC', 'Description'] %>% head


# Plot together with GM results
gonames = c("cytokine activity", "inflammatory response", "cytokine-mediated signaling pathway", "chemotaxis")
toplot = rbind(gm_gores_audu[gm_gores_audu$Description %in% gonames,],
		ju_gores_audu[ju_gores_audu$Description %in% gonames,])

toplot$log10FDR = -log10(toplot$qvalue)
toplot$Description = factor(toplot$Description, levels = rev(unique(toplot$Description)))

pdf('JUGM_GO_Enrichment_AU_UP_DU.pdf', width = 16)
ggscatter(toplot, x = 'log10FDR', y = 'Description', color = 'CellType', palette = c('blue', 'red'),
		size = 'OddsRatio', xlim = c(0, 25)) +
	ylab('') + xlab(expression('-log'[10]*'(FDR)')) +
	geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
	scale_size_continuous(range = c(6,10)) + 
	theme(text=element_text(size=25, face = 'bold'), legend.position = 'right')
dev.off()



gm_inf_gns = gm_gores_audu[gm_gores_audu$Description == 'inflammatory response', 'geneID'] %>% strsplit(., '/') %>% unlist
ju_inf_gns = ju_gores_audu[ju_gores_audu$Description == 'inflammatory response', 'geneID'] %>% strsplit(., '/') %>% unlist


setdiff(ju_inf_gns, gm_inf_gns) %>% as.data.frame %>% rio::export(., 'tmp.xlsx')
setdiff(gm_inf_gns, ju_inf_gns) %>% as.data.frame %>% rio::export(., 'tmp2.xlsx')



# Gene Ontology Enrichment AU-DU DOWN
au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_down = au_du[au_du$log2FoldChange < -0.6 & au_du$padj < 0.05,] %>% rownames
bcg = rownames(au_du)

gores_audu = GOenrich(au_du_down, bcg, qCut = 0.05)
gores_audu$ObsOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$ObsAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'GeneRatio']) %>% as.numeric()})
gores_audu$BcgOv = sapply(1:nrow(gores_audu), function(x){gsub('/.*', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$BcgAll = sapply(1:nrow(gores_audu), function(x){gsub('.*/', '', gores_audu[x, 'BgRatio']) %>% as.numeric()})
gores_audu$OddsRatio = (gores_audu$ObsOv / gores_audu$ObsAll) / (gores_audu$BcgOv / gores_audu$BcgAll)

rio::export(gores_audu, "~/project/Gozde_data/RNA-seq/files/JU_GO_Enrichment_AU_LOWER_DU.xlsx")
ju_gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/JU_GO_Enrichment_AU_LOWER_DU.xlsx")
gm_gores_audu = rio::import("~/project/Gozde_data/RNA-seq/files/GM_GO_Enrichment_AU_LOWER_DU.xlsx")

ju_gores_audu$CellType = 'Jurkat'
gm_gores_audu$CellType = 'GM'

gores_audu[gores_audu$GO == 'MF', 'Description'] %>% head
gores_audu[gores_audu$GO == 'BP', 'Description'] %>% head
gores_audu[gores_audu$GO == 'CC', 'Description'] %>% head

gonames = c("ribonucleoprotein complex binding", "ribosome biogenesis", "rRNA processing", "ribonucleoprotein complex biogenesis")
toplot = rbind(gm_gores_audu[gm_gores_audu$Description %in% gonames,],
		ju_gores_audu[ju_gores_audu$Description %in% gonames,])
toplot$log10FDR = -log10(toplot$qvalue)
toplot$Description = factor(toplot$Description, levels = rev(unique(toplot$Description)))

pdf('JUGM_GO_Enrichment_AU_DOWN_DU.pdf', width = 16)
ggscatter(toplot, x = 'log10FDR', y = 'Description', color = 'CellType', palette = c('blue', 'red'),
		size = 'OddsRatio', xlim = c(0, 40)) +
	ylab('') + xlab(expression('-log'[10]*'(FDR)')) +
	geom_vline(xintercept = 1.3, linetype = 'dashed', color = 'red') +
	scale_size_continuous(range = c(6,10)) + 
	theme(text=element_text(size=25, face = 'bold'), legend.position = 'right')
dev.off()



## HEATMAPS ##
load("~/project/Gozde_data/RNA-seq/R_objs/JU_counts_keep.RData")

JU_counts_norm = log2(cpm(JU_counts_keep) + 1)
JU_counts_norm_noai <- JU_counts_norm[, !(grepl(c("_AI"), colnames(JU_counts_norm)))]

JU_counts_scaled_noai <- scale(t(JU_counts_norm_noai))
JU_counts_scaled_noai <- t(JU_counts_scaled_noai)

# Significantly upregulated genes
au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames

au_di = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DI"), format = "DataFrame")
au_di = au_di[!(is.na(au_di$padj)),]
au_di_up = au_di[au_di$log2FoldChange > 0.6 & au_di$padj < 0.05,] %>% rownames

# Chemokine ligand genes #
chemogns = read.csv('~/project/Gozde_data/RNA-seq/files/Chemokine_Ligands_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

gns_toplot <- intersect(chemogns, rownames(JU_counts_scaled_noai))
#gns_toplot <- intersect(intersect(chemogns, au_du_up), au_di_up)

mat_plotgns <- JU_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3), paste0('JU_DI_', 1:3),
										paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_DUandAU_CCL_CXCLs_ALL_CCL1SIGN.pdf", width = 9)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()

# IL genes #

# Interleukin gene symbols in human
ilgns = read.csv('~/project/Gozde_data/RNA-seq/files/Interleukins_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

#gns_toplot <- intersect(ilgns, rownames(JU_counts_scaled_noai))
gns_toplot <- intersect(intersect(ilgns, au_du_up), au_di_up)

mat_plotgns <- JU_counts_scaled_noai[gns_toplot,]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3), paste0('JU_DI_', 1:3),
										paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_DUandAU_Interleukins_SIGN.pdf", width = 9)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()



# TNF genes #

# Interleukin gene symbols in human
tnfgns = read.csv('~/project/Gozde_data/RNA-seq/files/TumorNecrosisFactors_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

#gns_toplot <- intersect(tnfgns, rownames(JU_counts_scaled_noai))
gns_toplot <- intersect(intersect(tnfgns, au_du_up), au_di_up)
mat_plotgns <- JU_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3), paste0('JU_DI_', 1:3),
										paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_DUandAU_TNFs_SIGN.pdf", width = 9, height = 5)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()


# Cell Cycle Genes #
JU_counts_norm = log2(cpm(JU_counts_keep) + 1)
JU_counts_norm_audu <- JU_counts_norm[, !(grepl(c("_AI|_DI"), colnames(JU_counts_norm)))]

JU_counts_scaled_audu <- scale(t(JU_counts_norm_audu))
JU_counts_scaled_audu <- t(JU_counts_scaled_audu)

s_gns = intersect(cc.genes$s.genes, rownames(JU_counts_scaled_audu))
g2m_gns = intersect(cc.genes$g2m.genes, rownames(JU_counts_scaled_audu))

au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
s_gns_up = au_du[rownames(au_du) %in% s_gns & au_du$padj < 0.05 & au_du$log2FoldChange > 0, ] %>% rownames
s_gns_down = au_du[rownames(au_du) %in% s_gns & au_du$padj < 0.05 & au_du$log2FoldChange < 0, ] %>% rownames

g2m_gns_up = au_du[rownames(au_du) %in% g2m_gns & au_du$padj < 0.05 & au_du$log2FoldChange > 0, ] %>% rownames
g2m_gns_down = au_du[rownames(au_du) %in% g2m_gns & au_du$padj < 0.05 & au_du$log2FoldChange < 0, ] %>% rownames

# S
mat_plotgns <- JU_counts_scaled_audu[c(s_gns_up, s_gns_down),]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3), paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_S_Phase_Sign.pdf", width = 8, height = 10)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()

mat_plotgns <- JU_counts_norm_audu[c('CDC45', 'CDC6', 'PCNA', 'E2F8'),]
melted_int <- data.frame(melt(mat_plotgns))
melted_int$Group = gsub('_[0-9]', '', melted_int$Var2) %>% gsub('JU_', '', .)
melted_int$Group = factor(melted_int$Group, levels = c('DU', 'AU'))

pdf("Boxplot_JU_S_Phase_Sign_BoxPlot.pdf", width = 6, height = 8)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="Normalized Gene Expression")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  facet_wrap(~Var1, scales = 'free') + NoLegend() +
  rotate_x_text(45)
dev.off()


# G2M
mat_plotgns <- JU_counts_scaled_audu[c(g2m_gns_up, g2m_gns_down),]
melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3), paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_G2M_Phase_Sign.pdf", width = 8, height = 13)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = 0, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()



# Motif Enrichment Genes - Only AU and DU #
JU_counts_norm = log2(cpm(JU_counts_keep) + 1)
JU_counts_norm_noai <- JU_counts_norm[, !(grepl(c("_AI|_DI"), colnames(JU_counts_norm)))]
JU_counts_scaled_noai <- scale(t(JU_counts_norm_noai))
JU_counts_scaled_noai <- t(JU_counts_scaled_noai)

nfkb = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/nfkb_subunit_genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
fos = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/FOS_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
jun = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/JUN_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
bzip = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/bZIP_domain_genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']
irfs = read.csv('~/project/Gozde_data/ATACseq/files/Gene_Lists/IRF_Genes_HGNC.csv', skip = 1) %>% .[, 'Approved.symbol']

gns_toplot = c(nfkb, fos, jun, bzip, irfs)
gns_toplot <- intersect(gns_toplot, rownames(JU_counts_scaled_noai))

mat_plotgns <- JU_counts_scaled_noai[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
colnames(melted_int) <- c("gene", "sample", "norm_counts")

melted_int_ordered <- melted_int[order(melted_int$sample),]
melted_int_ordered$sample = factor(melted_int_ordered$sample, levels = c(paste0('JU_DU_', 1:3),
										paste0('JU_AU_', 1:3)))

pdf("heatMap_JU_DUandAU_MOTIFGENES_ALL.pdf", width = 9, height = 15)
ggplot(melted_int_ordered, aes(x=sample,y=gene,fill=norm_counts))+
  geom_tile(colour="white",size=0.2)+
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  rotate_x_text(45)
dev.off()


# Boxplots of motif enrichment genes
au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames
ju_au_du_up = au_du_up

gm_au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
gm_au_du = gm_au_du[!(is.na(gm_au_du$log2FoldChange)),]
gm_au_du_up = gm_au_du[gm_au_du$log2FoldChange > 0.6 & gm_au_du$padj < 0.05,] %>% rownames

JU_counts_norm = log2(cpm(JU_counts_keep) + 1)

gns_toplot <- intersect(gns_toplot, rownames(JU_counts_scaled_noai))
# plot only the genes significantly upregulated in AUvsDU in both GM and JU
gns_toplot <- c("FOS", "JUN", "ATF3", "CEBPB", "JDP2", "MAFF", "IRF1", "IRF7", "IRF9")
mat_plotgns <- JU_counts_norm[gns_toplot,]

melted_int <- data.frame(melt(mat_plotgns))
melted_int$Group = gsub('_[0-9]', '', melted_int$Var2) %>% gsub('JU_', '', .)

melted_int$Group = factor(melted_int$Group, levels = c('DU', 'DI', 'AU', 'AI'))
# change the level names for the paper
levels(melted_int$Group) <- c('CU', 'CI', 'AU', 'AI')

pdf("BoxPlot_JU_MotifGenes.pdf", width = 9, height = 12)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold')) +
  facet_wrap(~Var1, scales = 'free') + NoLegend() +
  rotate_x_text(45)
dev.off()


# plot the common motifs between GM and JU
pdf("~/project/Gozde_data/RNA-seq/files/BoxPlot_MotifGenes_common_JU.pdf", width = 9, height = 12)
ggboxplot(melted_int, x = 'Group', y = 'value', color = 'Group', add = 'jitter') +
  labs(x="",y="")+
  scale_fill_gradient2(midpoint = -0.5, low = 'blue', high = 'red')+
  theme_classic()+
  theme(text = element_text(size=25, face = 'bold'),legend.position = "none") +
  facet_wrap(~Var1, scales = 'free') + 
  rotate_x_text(45)
dev.off()


################################################################
## AI vs DU


ai_du = results(JU_RNA_DDS, c("Group_JU", "JU_AI", "JU_DU"), format = "DataFrame")
ai_du = ai_du[!(is.na(ai_du$padj)),]
ai_du_up = ai_du[ai_du$log2FoldChange > 0.6 & ai_du$padj < 0.05,] %>% rownames

au_du = results(JU_RNA_DDS, c("Group_JU", "JU_AU", "JU_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames

di_du = results(JU_RNA_DDS, c("Group_JU", "JU_DI", "JU_DU"), format = "DataFrame")
di_du = di_du[!(is.na(di_du$padj)),]
di_du_up = di_du[di_du$log2FoldChange > 0.6 & di_du$padj < 0.05,] %>% rownames

sum(ai_du_up %in% setdiff(au_du_up, di_du_up)) / length(ai_du_up)
sum(ai_du_up %in% setdiff(di_du_up, au_du_up)) / length(ai_du_up)

qw1 = setdiff(ai_du_up, union(di_du_up, au_du_up))


ai_du = results(GM_RNA_DDS, c("Group_GM", "GM_AI", "GM_DU"), format = "DataFrame")
ai_du = ai_du[!(is.na(ai_du$padj)),]
ai_du_up = ai_du[ai_du$log2FoldChange > 0.6 & ai_du$padj < 0.05,] %>% rownames

au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du = au_du[!(is.na(au_du$padj)),]
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames

di_du = results(GM_RNA_DDS, c("Group_GM", "GM_DI", "GM_DU"), format = "DataFrame")
di_du = di_du[!(is.na(di_du$padj)),]
di_du_up = di_du[di_du$log2FoldChange > 0.6 & di_du$padj < 0.05,] %>% rownames

sum(ai_du_up %in% setdiff(au_du_up, di_du_up)) / length(ai_du_up)
sum(ai_du_up %in% setdiff(di_du_up, au_du_up)) / length(ai_du_up)

qw2 = setdiff(ai_du_up, union(di_du_up, au_du_up))


setdiff(ai_du_up, union(di_du_up, au_du_up)) %>% as.data.frame %>% rio::export(., 'tmp.xlsx')



























