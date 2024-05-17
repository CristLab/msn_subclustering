library(parallel)
library(dplyr)
library(Seurat)
library(SoupX)
library(DoubletFinder)
library(ggplot2)
library(data.table)
library(BiocParallel)
library(SingleCellExperiment)
library(scDblFinder)


input.list <- scan("/home/crist/nac_rn7/nac_rn7_scdf_input.txt", what = list("character","character","numeric","character","character","character"))

rds.path <- input.list[[1]]
remove.samples <- input.list[[2]]
samples.to.remove <- unlist(strsplit(input.list[[3]],","))
path.to.save <- input.list[[4]]
metadata.path <- input.list[[5]]
project <- input.list[[6]]


sct <- readRDS(rds.path)
sct

head(sct@meta.data)
sct@meta.data$rowname <- rownames(sct@meta.data)
sct@meta.data$original_order <- 1:nrow(sct@meta.data) 
head(sct@meta.data)

Idents(sct) <- "subject"

#Remove samples with high contamination or other issues
if (remove.samples){
	sct <- subset(sct, idents = samples.to.remove, invert = TRUE)
	sct
	table(sct@meta.data$subject)
}

DefaultAssay(sct) <- "RNA"

#Run scDblFinder
sce <- as.SingleCellExperiment(sct)
sce <- scDblFinder(sce, samples="subject", BPPARAM=MulticoreParam(16), includePCs=20)

table(sce$seurat_clusters, sce$scDblFinder.class)

head(colData(sce))

doublet_calls <- as.data.frame(colData(sce)[ , c("cell","scDblFinder.class")])
head(doublet_calls)
sct@meta.data <- merge(sct@meta.data, doublet_calls, by=c("rowname" = "cell"))
# sct@meta.data <- merge(sct@meta.data, doublet_calls,by='row.names',all=TRUE)

head(sct@meta.data)

sct@meta.data <- sct@meta.data[order(sct@meta.data$original_order), ]

head(sct@meta.data)

rownames(sct@meta.data) <- sct@meta.data$rowname

head(sct@meta.data)

table(sct@meta.data$seurat_clusters, sct@meta.data$scDblFinder.class)

sct


#Change idents to doublet or singlet and print UMAP
Idents(sct) <- sct$scDblFinder.class
DimPlot(object = sct,reduction = "umap")

#Change idents back and print number of nuclei in different groups
Idents(sct) <- "seurat_clusters"
table(Idents(sct),sct@meta.data$scDblFinder.class)
table(sct@meta.data$subject,sct@meta.data$scDblFinder.class)

VlnPlot(sct, features = c("nCount_RNA"), pt.size=0, group.by = "scDblFinder.class")

saveRDS(sct, file = paste0(path.to.save,project,"_scdf.rds"))

sessionInfo()




