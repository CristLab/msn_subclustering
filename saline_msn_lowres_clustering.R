library(future)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2)

#Load full clean Seurat object
oudfull <- readRDS("/home/crist/OUD_10x/oudfull_clustered_clean.rds")
oudfull


#*****CHANGE THIS TO CLUSTER NAMES*****


#Subset to only MSN clusters
msn <- subset(x = oudfull, idents=c("15","16","17","18","20","21","22","23","24","25","26"))
msn

#Subset to just saline controls
Idents(msn) <- 'treatment'
saline.msn <- subset(x = msn, idents = c("saline"))
saline.msn

#Update Seurat object because original analysis was performed on Seurat v3
saline.msn <- UpdateSeuratObject(saline.msn)

#Run standard pipeline for normalizing and scaling data, finding variable features, and running PCA
#Note: this data set has already had mitochondrial transcripts removed as part of the original analysis
all.genes <- rownames(saline.msn)
saline.msn <- NormalizeData(saline.msn, normalization.method = "LogNormalize", scale.factor = 10000)
saline.msn <- ScaleData(saline.msn, features = all.genes, vars.to.regress = c("nCount_RNA"))
saline.msn <- FindVariableFeatures(saline.msn, selection.method = "mvp", mean.cutoff=c(0.003,2))
saline.msn <- RunPCA(saline.msn, features = VariableFeatures(object = saline.msn), npcs = 100)


#Low resolution clustering of MSNs
ElbowPlot(saline.msn, ndims = 100)

saline.msn <- FindNeighbors(saline.msn, dims = 1:20)
saline.msn <- FindClusters(saline.msn, resolution = 0.1)
saline.msn <- BuildClusterTree(saline.msn, reorder = TRUE, reorder.numeric = TRUE, dims = 1:20)

Idents(saline.msn) <- 'tree.ident'
levels(saline.msn) <- c(1,2,3,4,5,6)

PlotClusterTree(saline.msn)

#Print violin plot. QC check for clusters with abnormally high or low UMI distributions
VlnPlot(saline.msn, features = c("nCount_RNA"), pt.size=0, group.by='tree.ident') + theme(legend.position = 'none')

#Print table with number of nuclei from each subject in each cluster. QC check for biased distribution
table(Idents(saline.msn),saline.msn@meta.data$subject)

#Remove cluster 1 for abnormally low UMI number
saline.msn <- subset(x = saline.msn, idents=c("1"), invert = TRUE)

#Run UMAP
saline.msn <- RunUMAP(saline.msn,dims = 1:50, min.dist=0.5, n.epochs = 500)

p <- DimPlot(saline.msn, reduction = "umap", group.by='ident', shuffle=TRUE) + NoLegend()
p <- LabelClusters(plot = p, id = "ident", color = "black", repel=TRUE, box.padding = 2, segment.color = "black")
p <- p + theme(
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		axis.title=element_blank()
	)
p


saveRDS(saline.msn, file = "/home/crist/OUD_10x/Final_Objects/saline_msn_lowres_clean.rds")

