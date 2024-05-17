library(future)
library(dplyr)
library(Seurat)
library(ggplot2)

set.seed(123)

nac_rn7 <- readRDS("/home/crist/nac_rn7/nac_rn7_scdf.rds")
nac_rn7

#Keeping only singlets
Idents(nac_rn7) <- 'scDblFinder.class'
nac_rn7 <- subset(x = nac_rn7, idents=c("singlet"))
nac_rn7

#Cluster 7 has low UMI.
#Clusters 10, 15, 17, and 20 marked as majority doublet by scDblFinder
#Cluster 24 also has mixed markers by visual inspection
Idents(nac_rn7) <- 'seurat_clusters'
nac_rn7 <- subset(x = nac_rn7, idents=c("7","10","15","17","20","24"), invert = TRUE)
nac_rn7

#Make list of individual Seurat objects for each sample
seurat.list <- SplitObject(nac_rn7, split.by = "subject")

#SCTransform to get variable features and normalize individual objects
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


sct <- RunPCA(sct, npcs = 20, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:20, min.dist=0.5, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:20)
sct <- FindClusters(sct, resolution = 0.5)
sct

table(Idents(sct),sct@meta.data$subject)

#Violin plot of UMI for each cluster
VlnPlot(sct, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

#Plot starting UMAP
DimPlot(sct, reduction = "umap", group.by='ident')

#Dot plot for major markers using RNA assay
DefaultAssay(sct) <- "RNA"

sct <- NormalizeData(sct, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sct)
sct <- ScaleData(sct, features=all.genes, vars.to.regress = c("nCount_RNA"))

dot.features <- c("Vim","Cfap52","Gad2","Gad1","Slc17a8","Slc17a7","Slc17a6","Slc5a7","Th","Mog","Mag","Olig1","Olig2","Pcdh15","Pdgfra","Gfap","Ndrg2","Gja1","Aqp4","Pecam1","Cx3cr1","Mrc1","Drd3","Drd1","Drd2","Bcl11b","Pde1b")

p <- DotPlot(sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_preclean_dotplot1.png")

dot.features <- c("Gfral","Ret","Gdf15","Gcg","Glp1r","Gipr","Pirt","Pomc","Calcr","Calcrl","Ramp1","Ramp2","Ramp3","Lepr","Trpv1","Npy1r","Npy2r","Npy5r","Vip","Vipr1","Vipr2")

p <- DotPlot(sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_preclean_dotplot2.png")



###Remove residual problematic clusters. Cluster 20 for mixed markers and 22 for low UMI
Idents(sct) <- 'seurat_clusters'
sct <- subset(x = sct, idents=c("20","22"), invert = TRUE)
sct

table(Idents(sct),sct@meta.data$subject)

#Violin plot of UMI for each cluster
VlnPlot(sct, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

#Plot starting UMAP
DimPlot(sct, reduction = "umap", group.by='ident')

dot.features <- c("Vim","Cfap52","Gad2","Gad1","Slc17a8","Slc17a7","Slc17a6","Slc5a7","Th","Mog","Mag","Olig1","Olig2","Pcdh15","Pdgfra","Gfap","Ndrg2","Gja1","Aqp4","Pecam1","Cx3cr1","Mrc1","Drd3","Drd1","Drd2","Bcl11b","Pde1b")

p <- DotPlot(sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_clean_dotplot1.png")

dot.features <- c("Gfral","Ret","Gdf15","Gcg","Glp1r","Gipr","Pirt","Pomc","Calcr","Calcrl","Ramp1","Ramp2","Ramp3","Lepr","Trpv1","Npy1r","Npy2r","Npy5r","Vip","Vipr1","Vipr2")

p <- DotPlot(sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_clean_dotplot2.png")


VlnPlot(sct, features = c("Bcl11b"), pt.size=0) + theme(legend.position = 'none')

###Remove cluster 19 due to abnormally low UMI for neurons and mixed D1/D2 markers
sct <- subset(x = sct, idents=c("19"), invert=TRUE)
sct

saveRDS(sct, file = "/home/crist/nac_rn7/nac_rn7_clean.rds")

sessionInfo()