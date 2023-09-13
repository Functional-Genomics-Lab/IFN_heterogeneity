rm(list = ls())
library(dplyr)
library(ggpubr)
library(Seurat)
library(reshape2)
library(destiny)
library(ggrastr)

# custom function for the stacked bar plots
stackedbarplot = function(meta, groupx, groupfill, fn, horizontal = F, bold = F, wd = 10, hg = 10, adall = T){

require(tidyverse)
require(dplyr)
require(ggplot2)

clkeys = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_keys %>% as.data.frame
clkeys$size = meta %>% group_by(meta[, c(groupx, groupfill)]) %>% group_size 

# Add all cell numbers
if(adall == T){
	spkeys = meta %>% group_by(meta[,groupfill]) %>% group_keys %>% as.data.frame
	colnames(spkeys) = groupfill
	spkeys$size = meta %>% group_by(meta[, groupfill]) %>% group_size
	spkeys = spkeys %>% add_column(tmp = 'AllCells')
	colnames(spkeys)[3] = groupx
	clkeys = rbind(clkeys, spkeys)
}

colnames(clkeys) = c('cluster', 'variable', 'value')

plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
	    geom_bar(position = "fill", stat = "identity") +
	    xlab("") +
	    ylab("Percentage") +
	    theme_classic() +
	    #scale_fill_manual(values = colors) +
	    scale_y_continuous(labels = scales::percent_format()) +
	    theme(text=element_text(size=30), axis.text.y = element_text(face = 'bold')) +
	    rotate_x_text(45) + coord_flip()

if(horizontal == T){
plt = ggplot(clkeys, aes(x = cluster, y = value, fill = variable)) + 
	    geom_bar(position = "fill", stat = "identity") +
	    xlab("") +
	    ylab("Percentage") +
	    theme_classic() +
	    #scale_fill_manual(values = colors) +
	    scale_y_continuous(labels = scales::percent_format()) +
	    theme(text=element_text(size=30)) +
	    rotate_x_text(45)
}

pdf(paste0(fn, '.pdf'), width = wd, height = hg)
print(plt)
dev.off()

	return(plt)

}

# Load dataset
load('~/project/Gozde_data/scRNA-seq/in-house/Robjs/sub_seurat.RData')

dmso_seurat = subset(sub_seurat, subset = orig.ident == 'GM_SeV_DMSO')
dmso_seurat = NormalizeData(dmso_seurat)

# Analyze DMSO only #
dmso_seurat2 = dmso_seurat

# Infection percentage
dmso_seurat2[["inf_percent"]] = PercentageFeatureSet(object = dmso_seurat2, features = c('P', 'N', 'M', 'L', 'HN', 'F'))
#dmso_seurat2[["inf_percent_scaled"]] = scale(dmso_seurat2[["inf_percent"]]) 

#dmso_seurat2[["inf_status_2"]] = ifelse(dmso_seurat2[["inf_percent"]] > 0, 'Inf', 'Not_Inf')

#pdf("DMSO_Only_Inf_Status_2.pdf")
#DimPlot(dmso_seurat2, group.by = 'inf_status_2', pt.size = 1, raster = T) +
#theme(text=element_text(size=20, face = 'bold'), legend.position = "top") + ggtitle('')
#dev.off()




# Dimensionality reduction and clustering
dmso_seurat2 = SCTransform(dmso_seurat2, ncells = 1000)
dmso_seurat2 = RunPCA(dmso_seurat2, verbose = FALSE)
dmso_seurat2 = RunUMAP(dmso_seurat2, dims = 1:15)
dmso_seurat2 = FindNeighbors(dmso_seurat2, dims = 1:15, verbose = FALSE)
dmso_seurat2 = FindClusters(dmso_seurat2, verbose = FALSE, resolution = 0.1)

pdf("DMSO_Only_Inf_Percent.pdf")
FeaturePlot(dmso_seurat2, features = 'inf_percent', pt.size = 2, raster = T, sort = T) +
theme(text=element_text(size=20, face = 'bold'), legend.position = "top") + ggtitle('')
dev.off()

pdf("DMSO_Only_Inf_Status.pdf")
DimPlot(dmso_seurat2, group.by = 'inf_status', pt.size = 2, raster = T, cols = c('orange', 'blue')) +
theme(text=element_text(size=20, face = 'bold'), legend.position = "top") + ggtitle('')
dev.off()

pdf("DMSO_Only_IFN.pdf", width = 14)
FeaturePlot(dmso_seurat2, features = c('IFNB1', 'IFNL1'), sort = T, pt.size = 2, raster = T, keep.scale = 'all') +
theme(text=element_text(size=20, face = 'bold'))
dev.off()

# Plot IFN genes by cell cycle
IFN_genes = intersect(IFN_genes, rownames(dmso_seurat2))
dmso_seurat2$Phase = factor(dmso_seurat2$Phase, levels = c('G1', 'S', 'G2M'))

pdf('Broad_CellCycle_All_IFN.pdf', width = 8, height = 6)
DotPlot(dmso_seurat2, features = IFN_genes, group.by = 'Phase', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()

pdf('Broad_CellCycle_IFNB1_IFNL1.pdf', width = 8, height = 4)
DotPlot(dmso_seurat2, features = c('IFNB1', 'IFNL1'), group.by = 'Phase', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()




# Cell cycle genes from seurat
cc_genes = c(cc.genes$s.genes, cc.genes$g2m.genes)
cc_genes = intersect(cc_genes, rownames(dmso_seurat))

## Run Destiny ##
mat = dmso_seurat@assays$RNA@data[cc_genes, ]
#keepcells = sample(1:ncol(mat), 1000)
#mat2 = mat[, keepcells]
mat2 = mat
dm = DiffusionMap(t(as.matrix(mat2)), n_pcs = 50, k = 1000)

saveRDS(dm, '~/project/Gozde_data/scRNA-seq/in-house/Robjs/DMSO_Seurat_DiffusionMap.RDS')
dm = readRDS('~/project/Gozde_data/scRNA-seq/in-house/Robjs/DMSO_Seurat_DiffusionMap.RDS')

cmncells = intersect(names(dm$DC1), colnames(mat2))

dpm_df = data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  color = dmso_seurat$Phase[cmncells],
		  Gene1 = colSums( dmso_seurat@assays$RNA@data[c('IFNL1', 'IFNB1'), cmncells]) )

dpm_df$color = factor(dpm_df$color, levels = c('G1', 'S', 'G2M'))

ggplt = ggplot(dpm_df, aes(x = DC1, y = DC2, color = color)) +
    scale_color_manual(values = c('blue', 'orange', 'red')) +
    geom_point() + #scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic() +
    theme(text=element_text(size=20, face = 'bold'), legend.pos = 'top')

pdf('Pseudotime_CellCycle_DMSO.pdf', width = 8)
rasterize(ggplt, layers='Point', dpi=300)
dev.off()


# Find the cell in the middle
midpoint = sum(max(dpm_df$DC1), min(dpm_df$DC1))/2

incdf = dpm_df[dpm_df$DC1 > midpoint,]
decdf = dpm_df[dpm_df$DC1 < midpoint,]
toG2M = incdf[order(incdf$DC1)[1:1000],]
toG1 = decdf[order(decdf$DC1, decreasing = T)[1:5000],]

# Plot clusters with borders
toplot = melt(dpm_df[,1:3], id.vars = c('DC1', 'DC2'))
toplot$value = factor(toplot$value, levels = c('G1', 'S', 'G2M'))

ggplt = ggscatter(toplot, x = 'DC1', y = 'DC2', color = 'value', palette = c('blue', 'orange', 'red')) +
xlab("Diffusion component 1") + 
ylab("Diffusion component 2") +
theme_classic() +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'top') +
geom_vline(xintercept = c(min(toG1$DC1), max(toG2M$DC1)), linetype = 'dashed', color = 'darkgreen')

pdf('Pseudotime_CellCycle_DMSO_GROUPED.pdf', width = 8)
rasterize(ggplt, layers='Point', dpi=300)
dev.off()


# These are the cells between the brackets
midCells = c(rownames(toG1), rownames(toG2M))
G2M = incdf[order(incdf$DC1)[1001:nrow(incdf)],] %>% rownames
G1 = decdf[order(decdf$DC1, decreasing = T)[5001:nrow(decdf)],] %>% rownames

g1Cells = dmso_seurat@meta.data[dmso_seurat$Phase == 'G1',] %>% rownames
g2mCells = dmso_seurat@meta.data[dmso_seurat$Phase == 'G2M',] %>% rownames

g1early = setdiff(g1Cells, midCells)
g1late = intersect(g1Cells, midCells)

g2mlate = setdiff(g2mCells, midCells)
g2mearly = intersect(g2mCells, midCells)

meta = dmso_seurat@meta.data
meta$dmannot = 'S'
meta[g1early, 'dmannot'] = 'G1_Early'
meta[g1late, 'dmannot'] = 'G1_Late'
meta[g2mearly, 'dmannot'] = 'G2M_Early'
meta[g2mlate, 'dmannot'] = 'G2M_Late'
#meta[meta$dmannot == 'S', 'Phase'] %>% table
dmso_seurat@meta.data = meta
dmso_seurat$dmannot = factor(dmso_seurat$dmannot, levels = c('G1_Early', 'G1_Late', 'S', 'G2M_Early', 'G2M_Late'))

qw = rownames(dmso_seurat)
qw = qw[grepl('^NFKB', qw)]

pdf('CellCycle_CXCL.pdf', width = 10, height = 5)
DotPlot(dmso_seurat, features = qw, group.by = 'dmannot', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()


# Plot IFN genes by cell cycle
pdf('CellCycle_IFN.pdf', width = 10, height = 5)
DotPlot(dmso_seurat, features = c('IFNB1', 'IFNL1'), group.by = 'dmannot', dot.scale = 15) +
xlab('') + ylab('') +
theme(text=element_text(size=20, face = 'bold'), legend.pos = 'right') +
coord_flip()
dev.off()

stackedbarplot(dmso_seurat[[]], groupx = 'dmannot', groupfill = 'inf_status', 'CellCycle_by_infection_status')


# Replot the pseudotime with the new groups
dpm_df2 = data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  color = dmso_seurat$dmannot[cmncells] )

pdf('Pseudotime_FINER_CellCycle_DMSO.pdf', width = 10)
ggplot(dpm_df2, aes(x = DC1, y = DC2, color = color)) +
    scale_color_manual(values = c('lightblue', 'darkblue', 'orange', 'pink', 'red')) +
    geom_point(alpha = 0.2) + #scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic() +
    theme(text=element_text(size=20, face = 'bold'), legend.pos = 'top')
dev.off()



# Statistical test of IFN genes - Finer Phases
allgns = rownames(sub_seurat)
IFN_genes = allgns[grepl("IFN", allgns)]
(IFN_genes = IFN_genes[!grepl("R", IFN_genes)])

ifn_stats = FindMarkers(dmso_seurat, ident.1 = 'G1_Late', only.pos = T,
		features = IFN_genes, group.by = 'dmannot',
		min.pct = 0, logfc.threshold = 0)


# Plot table
library(gridExtra)
ifn_stats_grob = tableGrob(ifn_stats)

pdf('IFN_G1Late_vs_All.pdf', width = 10)
grid.arrange(ifn_stats_grob)
dev.off()



# Statistical test of IFN genes - Broad Phases
allgns = rownames(sub_seurat)
IFN_genes = allgns[grepl("IFN", allgns)]
(IFN_genes = IFN_genes[!grepl("R", IFN_genes)])

ifn_stats = FindMarkers(dmso_seurat, ident.1 = 'G1', only.pos = T,
		features = IFN_genes, group.by = 'Phase',
		min.pct = 0, logfc.threshold = 0)


# Plot table
library(gridExtra)
ifn_stats_grob = tableGrob(ifn_stats)

pdf('IFN_G1_vs_All.pdf', width = 10)
grid.arrange(ifn_stats_grob)
dev.off()



all_stats = FindMarkers(dmso_seurat, ident.1 = 'G1_Late', only.pos = T,
		features = rownames(dmso_seurat), group.by = 'dmannot',
		min.pct = 0.02, logfc.threshold = 0.05)
saveRDS(all_stats, '~/project/Gozde_data/scRNA-seq/in-house/Robjs/Late_G1_Markers.RDS')

all_stats = cbind(Gene = rownames(all_stats), all_stats)
rio::export(all_stats, '~/project/Gozde_data/scRNA-seq/in-house/Robjs/Late_G1_Markers.xlsx')


all_stats_sign = all_stats[all_stats$p_val_adj < 0.05,]


DotPlot(dmso_seurat, features = c('MCM2', 'MCM4'), group.by = 'dmannot')


library(dplyr)
library(ggpubr)
library(Seurat)
library(reshape2)
library(destiny)
library(ggrastr)
library(DESeq2)
library(edgeR)
source("~/onlybiohpc/pr3/OUR_DATA/utility_functions.R")

load(file="/home2/s422159/project/Gozde_data/RNA-seq/R_objs/GM_RNA_DDS.RData")

au_du = results(GM_RNA_DDS, c("Group_GM", "GM_AU", "GM_DU"), format = "DataFrame")
au_du_up = au_du[au_du$log2FoldChange > 0.6 & au_du$padj < 0.05,] %>% rownames


qw = au_du[intersect(rownames(au_du), rownames(all_stats_sign)),]
qw2 = qw[qw$padj < 0.05 & qw$log2FoldChange > 0.6,]
qw2[order(qw2$padj),] %>% head

qw[qw$padj < 0.05 & qw$log2FoldChange > 0.6,] %>% rownames













