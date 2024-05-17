library(future)
library(dplyr)
library(Seurat)
library(ggplot2)

set.seed(123)

nac_rn7 <- readRDS("/home/crist/nac_rn7/nac_rn7_clean.rds")
nac_rn7

### Remove cluster 19 due to poor markers and low UMI
# nac_rn7 <- subset(x = nac_rn7, idents=c("19"), invert=TRUE)
# nac_rn7

saveRDS(nac_rn7, file = "/home/crist/nac_rn7/nac_rn7_clean.rds")

###Subset to the MSN clusters based on high Bcl11b and Pde1b expression
msn <- subset(x = nac_rn7, idents=c("0","2","5","6","7","8","12","15","16"))
msn

###Make list of individual Seurat objects for each sample
seurat.list <- SplitObject(msn, split.by = "subject")

###Get new PCs
plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2, future.seed=TRUE)

seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                            SCTransform(x, vst.flavor = "v2", vars.to.regress = c("nCount_RNA"), verbose = FALSE) %>%
							RunPCA(npcs = 20, verbose = FALSE)
                         })
						 
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)


####Check references are correct each time if used
anchors <- FindIntegrationAnchors(object.list = seurat.list,
									reference = c(1,2,6,7),
									normalization.method = "SCT",
     								anchor.features = features)

sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

###Run 100 PCs here so they can be used downstream when subclustering
sct <- RunPCA(sct, npcs = 100, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:20, min.dist=0.5, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:20)
sct <- FindClusters(sct, resolution = 0.05)
sct

table(Idents(sct),sct@meta.data$subject)

###Violin plot of UMI for each cluster
VlnPlot(sct, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

###Plot starting UMAP and feature plots for D1 and D2 markers
DimPlot(sct, reduction = "umap", group.by='ident')

DefaultAssay(sct) <- "RNA"
FeaturePlot(sct, features = "Drd1")
FeaturePlot(sct, features = "Ebf1")
FeaturePlot(sct, features = "Htr4")
FeaturePlot(sct, features = "Isl1")
FeaturePlot(sct, features = "Drd2")
FeaturePlot(sct, features = "Scube1")
FeaturePlot(sct, features = "Stk32a")
FeaturePlot(sct, features = "Ppm1e")

saveRDS(sct, file = "/home/crist/nac_rn7/nac_rn7_msn.rds")


###Subcluster the major clusters
cluster0 <- subset(x = sct, idents=c("0"))
DefaultAssay(cluster0) <- "integrated"
cluster0

cluster0 <- RunPCA(cluster0, npcs = 20, verbose = FALSE)
cluster0 <- RunUMAP(cluster0, reduction = "pca", dims = 1:20, min.dist=0.5, n.epochs=500, verbose = FALSE)
cluster0 <- FindNeighbors(cluster0, reduction = "pca", dims = 1:20)
cluster0 <- FindClusters(cluster0, resolution = 0.65)
head(cluster0@meta.data)

VlnPlot(cluster0, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

DimPlot(cluster0, reduction = "umap", group.by='ident', raster=FALSE) + NoLegend()

metadata0 <- cluster0@meta.data[ , c("seurat_clusters"), drop = FALSE]
head(metadata0)
metadata0 <- metadata0 %>% 
       dplyr::rename("subclusters" = "seurat_clusters")

####as.numeric was adding 1 before the additional increment of 1. Seems to be caused by the column being factors instead of characters to start
metadata0$subclusters<-as.numeric(as.character(metadata0$subclusters)) + 1
head(metadata0)

metadata0$subclusters<-sub("^","D1_Ebf1_",metadata0$subclusters)
metadata0$subclusters<-as.character(metadata0$subclusters)
head(metadata0)



cluster1 <- subset(x = sct, idents=c("1"))
DefaultAssay(cluster1) <- "integrated"
cluster1

cluster1 <- RunPCA(cluster1, npcs = 20, verbose = FALSE)
cluster1 <- RunUMAP(cluster1, reduction = "pca", dims = 1:20, min.dist=0.5, n.epochs=500, verbose = FALSE)
cluster1 <- FindNeighbors(cluster1, reduction = "pca", dims = 1:20)
cluster1 <- FindClusters(cluster1, resolution = 0.75)
head(cluster1@meta.data)

VlnPlot(cluster1, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

DimPlot(cluster1, reduction = "umap", group.by='ident', raster=FALSE) + NoLegend()

metadata1 <- cluster1@meta.data[ , c("seurat_clusters"), drop = FALSE]
head(metadata1)
metadata1 <- metadata1 %>% 
       dplyr::rename("subclusters" = "seurat_clusters")

####as.numeric is adding 1 before the additional increment of 1. Seems to be caused by the column being factors instead of characters to start
metadata1$subclusters<-as.numeric(as.character(metadata1$subclusters)) + 1
head(metadata1)

metadata1$subclusters<-sub("^","D2_Stk32a_",metadata1$subclusters)
metadata1$subclusters<-as.character(metadata1$subclusters)
head(metadata1)



cluster2 <- subset(x = sct, idents=c("2"))
DefaultAssay(cluster2) <- "integrated"
cluster2

cluster2 <- RunPCA(cluster2, npcs = 20, verbose = FALSE)
cluster2 <- RunUMAP(cluster2, reduction = "pca", dims = 1:20, min.dist=0.5, n.epochs=500, verbose = FALSE)
cluster2 <- FindNeighbors(cluster2, reduction = "pca", dims = 1:20)
cluster2 <- FindClusters(cluster2, resolution = 0.15) 
head(cluster2@meta.data)

VlnPlot(cluster2, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

DimPlot(cluster2, reduction = "umap", group.by='ident', raster=FALSE) + NoLegend()

metadata2 <- cluster2@meta.data[ , c("seurat_clusters"), drop = FALSE]
head(metadata2)
metadata2 <- metadata2 %>% 
       dplyr::rename("subclusters" = "seurat_clusters")

####as.numeric is adding 1 before the additional increment of 1. Seems to be caused by the column being factors instead of characters to start
metadata2$subclusters<-as.numeric(as.character(metadata2$subclusters)) + 1
head(metadata2)

metadata2$subclusters<-sub("^","D1_Ppm1e_",metadata2$subclusters)
metadata2$subclusters<-as.character(metadata2$subclusters)
head(metadata2)



cluster3 <- subset(x = sct, idents=c("3"))
DefaultAssay(cluster3) <- "integrated"
cluster3

cluster3 <- RunPCA(cluster3, npcs = 20, verbose = FALSE)
cluster3 <- RunUMAP(cluster3, reduction = "pca", dims = 1:20, min.dist=0.5, n.epochs=500, verbose = FALSE)
cluster3 <- FindNeighbors(cluster3, reduction = "pca", dims = 1:20)
cluster3 <- FindClusters(cluster3, resolution = 0.25)
head(cluster3@meta.data)

VlnPlot(cluster3, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

DimPlot(cluster3, reduction = "umap", group.by='ident', raster=FALSE) + NoLegend()

metadata3 <- cluster3@meta.data[ , c("seurat_clusters"), drop = FALSE]
head(metadata3)
metadata3 <- metadata3 %>% 
       dplyr::rename("subclusters" = "seurat_clusters")

####as.numeric is adding 1 before the additional increment of 1. Seems to be caused by the column being factors instead of characters to start
metadata3$subclusters<-as.numeric(as.character(metadata3$subclusters)) + 1
head(metadata3)

metadata3$subclusters<-sub("^","D1_Htr4_",metadata3$subclusters)
metadata3$subclusters<-as.character(metadata3$subclusters)
head(metadata3)



cluster4 <- subset(x = sct, idents=c("4"))
DefaultAssay(cluster4) <- "integrated"
cluster4

cluster4 <- RunPCA(cluster4, npcs = 20, verbose = FALSE)
cluster4 <- RunUMAP(cluster4, reduction = "pca", dims = 1:20, min.dist=0.5, n.epochs=500, verbose = FALSE)
cluster4 <- FindNeighbors(cluster4, reduction = "pca", dims = 1:20)
cluster4 <- FindClusters(cluster4, resolution = 0.15) #Previously res=0.25
head(cluster4@meta.data)

VlnPlot(cluster4, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

DimPlot(cluster4, reduction = "umap", group.by='ident', raster=FALSE) + NoLegend()

metadata4 <- cluster4@meta.data[ , c("seurat_clusters"), drop = FALSE]
head(metadata4)
metadata4 <- metadata4 %>% 
       dplyr::rename("subclusters" = "seurat_clusters")

####as.numeric is adding 1 before the additional increment of 1. Seems to be caused by the column being factors instead of characters to start
metadata4$subclusters<-as.numeric(as.character(metadata4$subclusters)) + 1
head(metadata4)

metadata4$subclusters<-sub("^","D2_Scube1_",metadata4$subclusters)
metadata4$subclusters<-as.character(metadata4$subclusters)
head(metadata4)


#Use rbind to concatenate all of the metadata dataframes
metadata <- do.call("rbind", list(metadata0, metadata1, metadata2, metadata3, metadata4))

#####Merge reorders the columns which later breaks the Idents in Seurat. This setup makes a new column with the original order, then later sorts by it after the merge
sct@meta.data$original_order <- 1:nrow(sct@meta.data) 
head(sct@meta.data)

sct@meta.data <-merge(sct@meta.data,metadata,by='row.names',all=TRUE)
head(sct@meta.data)
levels(sct@meta.data$subclusters)

sct@meta.data <- sct@meta.data[order(sct@meta.data$original_order), ]
head(sct@meta.data)

sct@meta.data$original_order <- NULL
head(sct@meta.data)

print("Number of NAs")
length(is.na(sct@meta.data$subclusters))

rownames(sct@meta.data) <- sct@meta.data$Row.names
sct@meta.data$Row.names <- NULL
head(sct@meta.data)


# DefaultAssay(sct) <- "RNA"

# sct <- NormalizeData(sct, normalization.method = "LogNormalize", scale.factor = 10000)

# Idents(sct) <- 'subclusters'
# markers <-FindAllMarkers(sct, only.pos=TRUE)
# write.csv(markers, file= "/home/crist/nac_rn7/nac_rn7_subclusters_preclean.findallmarkers.csv")

DefaultAssay(sct) <- "integrated"

Idents(sct) <- 'subclusters'
sct <- RunPCA(sct, npcs = 100, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:100)
sct <- RunUMAP(sct,dims = 1:100, n.neighbors=400, repulsion.strength = 10, min.dist=1, n.epochs=1000, return.model = TRUE)

p <- DimPlot(sct, reduction = "umap", group.by='ident', shuffle=TRUE) + NoLegend()
p <- LabelClusters(plot = p, id = "ident", color = "black", segment.color = "black")
p <- p + theme(
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		axis.title=element_blank()
	)
p

table(sct@meta.data$subclusters,sct@meta.data$subject)



####Remove four clusters. Two of them have nothing but weak oligodendrocyte/ribosomal genes in FindMarkers. The other two have no markers except weak oxphos genes
#Presumed residual doublets and low quality droplets
#Rename the higher cluster names in Ebf1 and Stk32a in Idents and the subclusters metadata column
Idents(sct) <- 'subclusters'
sct <- subset(x = sct, idents=c(
"D1_Ebf1_11",
"D1_Ebf1_13",
"D2_Stk32a_11",
"D2_Stk32a_13"), invert=TRUE)

newident <- Idents(sct)
newident <-gsub("^D1_Ebf1_12$","D1_Ebf1_11",newident)
newident <-gsub("^D1_Ebf1_14$","D1_Ebf1_12",newident)
newident <-gsub("^D2_Stk32a_12$","D2_Stk32a_11",newident)
newident <-gsub("^D2_Stk32a_14$","D2_Stk32a_12",newident)

Idents(sct) <- newident
sct@meta.data$subclusters <- as.character(Idents(sct))

table(sct@meta.data$subclusters,sct@meta.data$subject)

#"Randomize" levels for better UMAP coloring
levels(sct) <- c(
"D1_Ebf1_5",
"D2_Stk32a_1",
"D2_Stk32a_11",
"D1_Ebf1_4",
"D2_Stk32a_12",
"D1_Htr4_1",
"D2_Stk32a_7",
"D1_Ebf1_8",
"D2_Stk32a_9",
"D1_Ebf1_7",
"D2_Stk32a_3",
"D1_Ppm1e_1",
"D1_Ebf1_2",
"D2_Stk32a_8",
"D1_Ebf1_10",
"D1_Ppm1e_3",
"D1_Htr4_2",
"D1_Ebf1_1",
"D2_Stk32a_5",
"D1_Ebf1_6",
"D2_Scube1_1",
"D2_Stk32a_2",
"D1_Ebf1_9",
"D1_Htr4_3",
"D2_Stk32a_6",
"D1_Ebf1_12",
"D1_Ppm1e_4",
"D2_Stk32a_10",
"D1_Ebf1_3",
"D1_Ppm1e_2",
"D2_Stk32a_4",
"D2_Scube1_2",
"D1_Htr4_4",
"D1_Ebf1_11"
)



sct <- RunPCA(sct, npcs = 100, verbose = FALSE)
sct <- FindNeighbors(sct, dims = 1:100)
sct <- RunUMAP(sct,dims = 1:100, n.neighbors=400, repulsion.strength = 10, min.dist=1, n.epochs=500, return.model = TRUE)

p <- DimPlot(sct, reduction = "umap", group.by='ident', shuffle=TRUE) + NoLegend()
p <- LabelClusters(plot = p, id = "ident", color = "black", segment.color = "black")
p <- p + theme(
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		axis.title=element_blank()
	)
p

DefaultAssay(sct) <- "RNA"

sct <- NormalizeData(sct, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features=all.genes, vars.to.regress = c("nCount_RNA"))

saveRDS(sct, file = "/home/crist/nac_rn7/nac_rn7_subcluster.rds")

sessionInfo()