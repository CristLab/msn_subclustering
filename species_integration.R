library(future)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(nichenetr)

human.nac <- readRDS("/home/crist/nac_rn7/human_nac.rds")
human.nac

human.msn <- subset(x = human.nac, idents = c("MSN.D1_A","MSN.D1_B","MSN.D1_C","MSN.D1_D","MSN.D1_E","MSN.D1_F","MSN.D2_A","MSN.D2_B","MSN.D2_C","MSN.D2_D"))
human.msn
human.msn <- AddMetaData(human.msn, metadata = "human", col.name = 'species')
head(human.msn@meta.data)

mouse.msn <- readRDS("/home/crist/nac_rn7/mouse_msn_clean.rds")
mouse.msn
mouse.msn <- AddMetaData(mouse.msn, metadata = "mouse", col.name = 'species')
head(mouse.msn@meta.data)

rat.msn <- readRDS("/home/crist/nac_rn7/nac_rn7_subcluster.rds")
rat.msn
rat.msn <- AddMetaData(rat.msn, metadata = "rat", col.name = 'species')
head(rat.msn@meta.data)

DefaultAssay(rat.msn) <- "RNA"

#Code to convert from human to mouse genes is here: https://github.com/satijalab/seurat/issues/2617

exp_mtx <- as.matrix(human.msn@assays$originalexp@data)

con_df <- data.frame(hum_orig = rownames(exp_mtx),
                     mouse = convert_human_to_mouse_symbols(rownames(exp_mtx)),
                     stringsAsFactors = F)

## Remove NAs where there is no match
con_df <- con_df[!is.na(con_df$mouse),,F]

## Filter the expression matrix for genes which a mouse counterpart is available
exp_mtx <- exp_mtx[con_df$hum_orig,]

## Now change the rownames of the matrix to the mouse gene names
rownames(exp_mtx) <- con_df$mouse

## Create the seurat object with mouse genes.
human.msn.mouse <- CreateSeuratObject(counts = exp_mtx, meta.data = human.msn@meta.data )



plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2)

human.msn.mouse <- NormalizeData(human.msn.mouse, normalization.method = "LogNormalize", scale.factor = 10000)
human.msn.mouse <- FindVariableFeatures(human.msn.mouse, selection.method = "mvp", mean.cutoff=c(0.003,2))
human.msn.mouse

total.genes <- list(rownames(human.msn.mouse),
                    rownames(rat.msn),
                    rownames(mouse.msn))

common.genes <- Reduce(f = intersect, x = total.genes)

rat.msn <- ScaleData(rat.msn, verbose = FALSE)
rat.msn <- RunPCA(rat.msn, features = common.genes, verbose = FALSE)

mouse.msn <- ScaleData(mouse.msn, verbose = FALSE)
mouse.msn <- RunPCA(mouse.msn, features = common.genes, verbose = FALSE)

human.msn.mouse <- ScaleData(human.msn.mouse, verbose = FALSE)
human.msn.mouse <- RunPCA(human.msn.mouse, features = common.genes, verbose = FALSE)


features <- SelectIntegrationFeatures(object.list = c(rat.msn, mouse.msn, human.msn.mouse))
###Using rpca instead of cca based on some forum comments about relative weighting of batch correction vs bio-conservation
anchors <- FindIntegrationAnchors(object.list = c(rat.msn, mouse.msn, human.msn.mouse), dims=1:10,anchor.features = features, reduction = "rpca")
msn.combined <- IntegrateData(anchorset = anchors, dims=1:10)
msn.combined

DefaultAssay(msn.combined) <- "integrated"

msn.combined <- RunPCA(msn.combined, npcs = 10, verbose = FALSE)
msn.combined <- RunUMAP(msn.combined, reduction = "pca", dims = 1:10, min.dist=0.5, n.epochs=500)

msn.list <- SplitObject(msn.combined, split.by = "species")

DimPlot(msn.combined, reduction = "umap", shuffle=TRUE, raster=FALSE, group.by="species") + NoAxes()

#UMAP coloring by original called clusters for each species
DimPlot(msn.list[["human"]], reduction = "umap", shuffle=TRUE, raster=FALSE, group.by="cellType") + NoAxes() + NoLegend()

DimPlot(msn.list[["rat"]], reduction = "umap", shuffle=TRUE, raster=FALSE, group.by="seurat_clusters") + NoAxes() + NoLegend()

DimPlot(msn.list[["mouse"]], reduction = "umap", shuffle=TRUE, raster=FALSE, group.by="orig.ident") + NoAxes() + NoLegend()

saveRDS(msn.combined, file = "/home/crist/nac_rn7/cross_species_msn.rds")

sessionInfo()

