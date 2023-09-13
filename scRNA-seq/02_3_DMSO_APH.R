rm(list = ls())
library(dplyr)
library(ggpubr)
library(Seurat)
library(reshape2)
library(destiny)
library(ggrastr)
library(gridExtra)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

# Load dataset
load('~/project/Gozde_data/scRNA-seq/in-house/Robjs/sub_seurat.RData')


####
## QC PLOTS
####

# Plot UMI and cell number
meta = sub_seurat@meta.data
meta$log10UMI = log10(meta$nCount_RNA)

pdf('DMSO_APH_Number_of_UMI.pdf')
ggboxplot(meta, x = 'orig.ident', y = 'log10UMI', color = 'orig.ident',
		palette = c('blue', 'red'), ylim = c(2,5)) +
    xlab("") + 
    ylab("Number of UMI (log10)") +
    theme_classic() +
    theme(text=element_text(size=20, face = 'bold'), legend.pos = 'top') +
    NoLegend()
dev.off()


dmso_cells = sum(meta$orig.ident == 'GM_SeV_DMSO')
aph_cells = sum(meta$orig.ident == 'GM_SeV_APH')
toplot = data.frame(vars = c('DMSO', 'APH'), vals = c(dmso_cells, aph_cells))
pdf('DMSO_APH_Number_of_Cells.pdf')
ggbarplot(toplot, x = 'vars', y = 'vals', color = 'vars', fill = 'vars',
		palette = c('blue', 'red'), ylim = c(0,12000)) +
    xlab("") + 
    ylab("Number of cells") +
    theme_classic() +
    theme(text=element_text(size=20, face = 'bold'), legend.pos = 'top') +
    NoLegend()
dev.off()

# UMAP without harmony
sub_seurat2 = sub_seurat

sub_seurat2 = SCTransform(sub_seurat2, ncells = 1000)
sub_seurat2 = RunPCA(sub_seurat2, verbose = FALSE)
sub_seurat2 = RunUMAP(sub_seurat2, dims = 1:15)

pdf("DMSO_APH_Before_Harmony.pdf")
DimPlot(sub_seurat2, group.by = 'orig.ident', cols = c('red', 'blue'), raster = T)
dev.off()

pdf("DMSO_APH_After_Harmony.pdf")
DimPlot(sub_seurat, group.by = 'orig.ident', cols = c('red', 'blue'), raster = T)
dev.off()



####
## GROUP CLUSTERS
####

# Set Idents as clusters
Idents(sub_seurat) <- sub_seurat$old.ident
gr1 <- subset(sub_seurat, subset = old.ident %in% c("0", "1", "2", "4", "8"))
gr2 <- subset(sub_seurat, subset = old.ident %in% c("3",  "9"))
gr3 <- subset(sub_seurat, subset = old.ident %in% c("5",  "6"))
gr4 <- subset(sub_seurat, subset = old.ident %in% c("7",  "10"))

# create a new meta column which indicates these custom groups
tmp <- ifelse(colnames(sub_seurat) %in% colnames(gr1), "Group1", ifelse(colnames(sub_seurat) %in% colnames(gr2), "Group2", ifelse(colnames(sub_seurat) %in% colnames(gr3), "Group3", "Group4")))

table(tmp)
sub_seurat$custom_gr <- tmp

# Colors
library(RColorBrewer)
cols = brewer.pal(n = 4, name = 'Dark2')

pdf('DMSO_APH_Custom_Groups.pdf')
DimPlot(sub_seurat, label = T, raster = T, group.by = 'custom_gr', label.size = 7, cols = cols) +
NoLegend() + ggtitle('')
dev.off()

# Infection status
stackedbarplot(sub_seurat[[]], groupx='custom_gr', groupfill='inf_status', fn = 'stackeBbarPlot_custom_gr_infected_cell_composition')

# DMSO/APH composition
stackedbarplot(sub_seurat[[]], groupx='custom_gr', groupfill='orig.ident', fn = 'stackeBbarPlot_custom_gr_APH_DMSO_composition')


pdf('AllGroups_IFNB1_IFNL1.pdf', width = 8, height = 5)
DotPlot(sub_seurat, features = c('IFNB1', 'IFNL1'), group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()


sub_seurat_nogr4 = subset(sub_seurat, subset = custom_gr == 'Group4', invert = T)
pdf('NoGroup4_IFNB1_IFNL1.pdf', width = 8, height = 5)
DotPlot(sub_seurat_nogr4, features = c('IFNB1', 'IFNL1'), group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()


# Plot viral mRNAs across groups
vir_gns = c("N",  "P",  "M",  "F",  "HN", "L" )

pdf('DMSO_APH_Viral_Genes_Dotplot.pdf', width = 10, height = 5)
DotPlot(sub_seurat, features = vir_gns, group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()


# Plot s genes across groups
s_gns = c("CDC6", "CDC45", "E2F8", "PCNA", "CCNE2")

pdf('DMSO_APH_S_Genes_Dotplot.pdf', width = 10, height = 5)
DotPlot(sub_seurat, features = s_gns, group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()




# check the expression of all IFN genes
allgns = rownames(sub_seurat)
IFN_genes <- allgns[grepl("IFN", allgns)]
(IFN_genes <- IFN_genes[!grepl("R", IFN_genes)])


pdf('AllGroups_AllIFN.pdf', width = 8, height = 5)
DotPlot(sub_seurat, features = IFN_genes, group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
rotate_x_text(90)
dev.off()

pdf('NoGroup4_IFNB1_IFNL1.pdf', width = 8, height = 5)
DotPlot(sub_seurat_nogr4, features = IFN_genes, group.by = 'custom_gr', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
rotate_x_text(90)
dev.off()


# Statistical test of IFN genes
allgns = rownames(sub_seurat)
IFN_genes = allgns[grepl("IFN", allgns)]
(IFN_genes = IFN_genes[!grepl("R", IFN_genes)])

gr2_gr1 = FindMarkers(sub_seurat, ident.1 = 'Group2', ident.2 = 'Group1', only.pos = T,
		features = IFN_genes, group.by = 'custom_gr',
		min.pct = 0, logfc.threshold = 0)

gr3_gr1 = FindMarkers(sub_seurat, ident.1 = 'Group3', ident.2 = 'Group1', only.pos = T,
		features = IFN_genes, group.by = 'custom_gr',
		min.pct = 0, logfc.threshold = 0)

gr4_gr1 = FindMarkers(sub_seurat, ident.1 = 'Group4', ident.2 = 'Group1', only.pos = T,
		features = IFN_genes, group.by = 'custom_gr',
		min.pct = 0, logfc.threshold = 0)


gr2_gr1_grob = tableGrob(gr2_gr1)
gr3_gr1_grob = tableGrob(gr3_gr1)
gr4_gr1_grob = tableGrob(gr4_gr1)

pdf('IFN_Group2_vs_Group1.pdf', width = 10)
grid.arrange(gr2_gr1_grob)
dev.off()

pdf('IFN_Group3_vs_Group1.pdf', width = 10)
grid.arrange(gr3_gr1_grob)
dev.off()

pdf('IFN_Group4_vs_Group1.pdf', width = 10)
grid.arrange(gr4_gr1_grob)
dev.off()

saveRDS(sub_seurat, '~/project/Gozde_data/scRNA-seq/in-house/Robjs/sub_seurat_grouped.RDS')




