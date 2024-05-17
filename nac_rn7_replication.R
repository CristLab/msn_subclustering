library(future)
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(clustifyr)
library(pheatmap)
library(RColorBrewer)
library(SCpubr)
library(patchwork)
library(ggplotify)

saline.msn <- readRDS("/home/crist/nac_rn7/nac_rn7_subcluster.rds")
saline.msn

DefaultAssay(saline.msn) <- "RNA"
saline.msn <- FindVariableFeatures(saline.msn, selection.method = "mvp", mean.cutoff=c(0.003,2))

#Read in the count matrices from the two publications from the Phillips et al and create Seurat objects
#Add percentage of mitochondrial transcripts to metadata of each object and remove all mitochondrial transcripts
#Generate percentile tables for number of features and UMI
#Remove nuclei with low number of genes, high number of UMI, and/or mitochondrial fraction 5% or more
srm.data <- Read10X(data.dir="/home/crist/OUD_10x/Day_NAc/Saline_Repeated_Male/")
srm <- CreateSeuratObject(counts = srm.data, min.cells=1)
rm(srm.data)
saveRDS(srm, file = "/home/crist/nac_rn7/salinerepeatedmale.rds")

srm[["percent.mt"]] <- PercentageFeatureSet(srm, pattern = "^Mt-")

all.genes <- as.data.frame(rownames(srm))
colnames(all.genes) <- c("genes")
nomtgenes <- all.genes[!grepl("^Mt-",all.genes$genes),]
srm <- subset(srm, features = nomtgenes)

nfeatures <- srm@meta.data$nFeature_RNA
quantile(nfeatures, probs = seq(0,1,0.005))

ncounts <- srm@meta.data$nCount_RNA
quantile(ncounts, probs = seq(0,1,0.005))

srm <- subset(srm, subset = nFeature_RNA > 499 & nCount_RNA < 19428 & percent.mt < 5)


srf.data <- Read10X(data.dir="/home/crist/OUD_10x/Day_NAc/Saline_Repeated_Female/")
srf <- CreateSeuratObject(counts = srf.data, min.cells=1)
rm(srf.data)
saveRDS(srf, file = "/home/crist/nac_rn7/salinerepeatedfemale.rds")

srf[["percent.mt"]] <- PercentageFeatureSet(srf, pattern = "^Mt-")

all.genes <- as.data.frame(rownames(srf))
colnames(all.genes) <- c("genes")
nomtgenes <- all.genes[!grepl("^Mt-",all.genes$genes),]
srf <- subset(srf, features = nomtgenes)

nfeatures <- srf@meta.data$nFeature_RNA
quantile(nfeatures, probs = seq(0,1,0.005))

ncounts <- srf@meta.data$nCount_RNA
quantile(ncounts, probs = seq(0,1,0.005))

srf <- subset(srf, subset = nFeature_RNA > 499 & nCount_RNA < 21844 & percent.mt < 5)


sam.data <- Read10X(data.dir="/home/crist/OUD_10x/Day_NAc/Saline_Acute_Male/")
sam <- CreateSeuratObject(counts = sam.data, min.cells=1)
rm(sam.data)
saveRDS(sam, file = "/home/crist/nac_rn7/salineacutemale.rds")

sam[["percent.mt"]] <- PercentageFeatureSet(sam, pattern = "^Mt-")

all.genes <- as.data.frame(rownames(sam))
colnames(all.genes) <- c("genes")
nomtgenes <- all.genes[!grepl("^Mt-",all.genes$genes),]
sam <- subset(sam, features = nomtgenes)

nfeatures <- sam@meta.data$nFeature_RNA
quantile(nfeatures, probs = seq(0,1,0.005))

ncounts <- sam@meta.data$nCount_RNA
quantile(ncounts, probs = seq(0,1,0.005))

sam <- subset(sam, subset = nFeature_RNA > 499 & nCount_RNA < 21279 & percent.mt < 5)


saf.data <- Read10X(data.dir="/home/crist/OUD_10x/Day_NAc/Saline_Acute_Female/")
saf <- CreateSeuratObject(counts = saf.data, min.cells=1)
rm(saf.data)
saveRDS(saf, file = "/home/crist/nac_rn7/salineacuteemale.rds")

saf[["percent.mt"]] <- PercentageFeatureSet(saf, pattern = "^Mt-")

all.genes <- as.data.frame(rownames(saf))
colnames(all.genes) <- c("genes")
nomtgenes <- all.genes[!grepl("^Mt-",all.genes$genes),]
saf <- subset(saf, features = nomtgenes)

nfeatures <- saf@meta.data$nFeature_RNA
quantile(nfeatures, probs = seq(0,1,0.005))

ncounts <- saf@meta.data$nCount_RNA
quantile(ncounts, probs = seq(0,1,0.005))

saf <- subset(saf, subset = nFeature_RNA > 499 & nCount_RNA < 23555 & percent.mt < 5)


#Merge samples into single Seurat object
replication.nac <- merge(srm, y = c(srf,sam,saf))
replication.nac

#Add metadata
subject <- sapply(strsplit(rownames(replication.nac@meta.data), split="_"), "[[", 2)
subject <-gsub("^1$","srm",subject)
subject <-gsub("^2$","srf",subject)
subject <-gsub("^3$","sam",subject)
subject <-gsub("^4$","saf",subject)
replication.nac <- AddMetaData(object = replication.nac, metadata=data.frame(subject=subject, row.names=rownames(replication.nac@meta.data)))

Idents(replication.nac) <- 'subject'

sex <- sapply(strsplit(rownames(replication.nac@meta.data), split="_"), "[[", 2)
sex <-gsub("^1$","male",sex)
sex <-gsub("^2$","female",sex)
sex <-gsub("^3$","male",sex)
sex <-gsub("^4$","female",sex)
replication.nac <- AddMetaData(object = replication.nac, metadata=data.frame(sex=sex, row.names=rownames(replication.nac@meta.data)))


group <- sapply(strsplit(rownames(replication.nac@meta.data), split="_"), "[[", 2)
group <-gsub("^1$","repeated",group)
group <-gsub("^2$","repeated",group)
group <-gsub("^3$","acute",group)
group <-gsub("^4$","acute",group)
replication.nac <- AddMetaData(object = replication.nac, metadata=data.frame(group=group, row.names=rownames(replication.nac@meta.data)))

head(replication.nac@meta.data)

#The acute controls are alligned to rn6 and the repeated controls are on rn7. 
#Subset to only genes that are in all samples to remove issues with gene symbol changes
total.genes <- list(rownames(srm),
                    rownames(srf),
                    rownames(sam),
                    rownames(saf))

common.genes <- Reduce(f = intersect, x = total.genes)
replication.nac <- subset(replication.nac, features = common.genes)
replication.nac

#Standard clustering pipeline
plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2)

replication.nac <- NormalizeData(replication.nac, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(replication.nac)
replication.nac <- ScaleData(replication.nac, features = all.genes)
replication.nac <- FindVariableFeatures(replication.nac, selection.method = "mvp", mean.cutoff=c(0.003,2))
replication.nac <- RunPCA(replication.nac, features = VariableFeatures(object = replication.nac))

#There are the parameters used in the Phillips et al 2023 manuscript
replication.nac <- FindNeighbors(replication.nac, dims = 1:17)
replication.nac <- FindClusters(replication.nac, resolution = 0.15)

#Dot plot for identifying MSNs
features <- c("Syt1","Gad1","Bcl11b","Foxp2","Rgs5","Cldn5","Arhgap15","Pdgfra","Gja1","Hapln2","Mbp","Slc17a7","Sst","Kit","Elavl2","Kcnc2","Grm8","Drd3","Drd2","Drd1","Ebf1","Slc5a7")

DotPlot(replication.nac, features=features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

#Violin plot for number of UMI. QC for removal of any MSN clusters that are low quality
VlnPlot(replication.nac, features = c("nCount_RNA"), pt.size=0, group.by='seurat_clusters') + theme(legend.position = 'none')

#Subset to just the MSN clusters. Based on high Bcl11b, high Foxp2, no Elavl2 expression
replication.msn <- subset(x = replication.nac, idents = c("0","2","5","12"))
replication.msn

replication.msn <- FindVariableFeatures(replication.msn, selection.method = "mvp", mean.cutoff=c(0.003,2))


DefaultAssay(saline.msn) <- "RNA"
saline.msn

saline.msn <- NormalizeData(saline.msn, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(saline.msn)
saline.msn <- ScaleData(saline.msn, features = all.genes)
saline.msn <- FindVariableFeatures(saline.msn, selection.method = "mvp", mean.cutoff=c(0.003,2))
saline.msn <- RunPCA(saline.msn, features = VariableFeatures(object = saline.msn), npcs = 100)

anchors <- FindTransferAnchors(reference = saline.msn, query = replication.msn,
    dims = 1:100, reference.reduction = "pca")

replication.msn <- MapQuery(anchorset = anchors, reference = saline.msn, query = replication.msn,
    refdata = list(celltype = "subclusters"), reference.dims=1:100, query.dims=1:100, reference.reduction = "pca", reduction.model = "umap")
	
saveRDS(replication.msn, file = "/home/crist/nac_rn7/replication_msn_clean.rds")

#"Randomize" the levels to minimize similar colors being assigned to adjacent clusters on the UMAP
Idents(saline.msn) <- 'subclusters'

levels(saline.msn) <- c(
"D1_Ebf1_5",
"D2_Stk32a_1",
"D2_Stk32a_11",
"D1_Ebf1_4",
"D2_Stk32a_12",
"D1_Htr4_1",
"D2_Stk32a_7",
"D1_Ebf1_8",
"D2_Stk32a_9",
"D1_Ebf1_12",
"D2_Stk32a_3",
"D1_Ppm1e_1",
"D1_Ebf1_2",
"D2_Stk32a_8",
"D1_Ebf1_10",
"D1_Ppm1e_3",
"D1_Htr4_2",
"D1_Ebf1_9",
"D2_Stk32a_5",
"D1_Ebf1_6",
"D2_Scube1_1",
"D2_Stk32a_2",
"D1_Ebf1_1",
"D1_Htr4_3",
"D2_Stk32a_6",
"D1_Ebf1_7",
"D1_Ppm1e_4",
"D2_Stk32a_10",
"D1_Ebf1_3",
"D1_Ppm1e_2",
"D2_Stk32a_4",
"D2_Scube1_2",
"D1_Htr4_4",
"D1_Ebf1_11"
)

templevels <- levels(saline.msn)

DefaultAssay(saline.msn) <- "RNA"
DefaultAssay(replication.msn) <- "RNA"


Idents(replication.msn) <- 'predicted.celltype'
levels(replication.msn) <- templevels

p1 <- DimPlot(saline.msn, reduction = "umap",shuffle = TRUE)+ ggtitle("Discovery") + theme_void() + NoLegend() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5)
		)
p2 <- DimPlot(replication.msn, reduction = "ref.umap",shuffle = TRUE) + ggtitle("Replication") + theme_void() + NoLegend() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5)
		)
	
p3 <- Nebulosa::plot_density(saline.msn, "Fermt1", reduction = "umap", method = "ks") + theme_void() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5),
		legend.key.width= unit(5, 'mm'),
		legend.title = element_text(size=8),
		legend.text = element_text(size=6)
		)

p4 <- Nebulosa::plot_density(saline.msn, "Col14a1", reduction = "umap", method = "ks") + theme_void() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5),
		legend.key.width= unit(5, 'mm'),
		legend.title = element_text(size=8),
		legend.text = element_text(size=6)
		)
	
p5 <- Nebulosa::plot_density(replication.msn, "Fermt1", reduction = "ref.umap", method = "ks") + theme_void() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5),
		legend.key.width= unit(5, 'mm'),
		legend.title = element_text(size=8),
		legend.text = element_text(size=6)
		)

p6 <- Nebulosa::plot_density(replication.msn, "Col14a1", reduction = "ref.umap", method = "ks") + theme_void() + theme(
		plot.title = element_text(hjust = 0.5,vjust = 0.5),
		legend.key.width= unit(5, 'mm'),
		legend.title = element_text(size=8),
		legend.text = element_text(size=6)
		)

#Use clustifyr to get correlation plot
s_ref <- seurat_ref(
  seurat_object = saline.msn,
  cluster_col = "subclusters"
)

replication.msn.meta <- replication.msn@meta.data

res <- clustify(
  input = replication.msn,
  metadata = replication.msn.meta,
  cluster_col = "predicted.celltype",
  ref_mat = s_ref,
  obj_out = FALSE,
  query_genes = VariableFeatures(saline.msn),
  compute_method = "pearson"
)

cluster_order <- c(
"D1_Ppm1e_1",
"D1_Ppm1e_2",
"D1_Ppm1e_3",
"D1_Ppm1e_4",
"D1_Htr4_1",
"D1_Htr4_2",
"D1_Htr4_3",
"D1_Htr4_4",
"D1_Ebf1_1",
"D1_Ebf1_2",
"D1_Ebf1_3",
"D1_Ebf1_4",
"D1_Ebf1_5",
"D1_Ebf1_6",
"D1_Ebf1_7",
"D1_Ebf1_8",
"D1_Ebf1_9",
"D1_Ebf1_10",
"D1_Ebf1_11",
"D1_Ebf1_12",
"D2_Scube1_1",
"D2_Scube1_2",
"D2_Stk32a_1",
"D2_Stk32a_2",
"D2_Stk32a_3",
"D2_Stk32a_4",
"D2_Stk32a_5",
"D2_Stk32a_6",
"D2_Stk32a_7",
"D2_Stk32a_8",
"D2_Stk32a_9",
"D2_Stk32a_10",
"D2_Stk32a_11",
"D2_Stk32a_12"
)

res <- res[order(match(row.names(res), cluster_order)),]
res <- res[,cluster_order]
res

cluster_group <- c(
"D1_Ppm1e",
"D1_Ppm1e",
"D1_Ppm1e",
"D1_Ppm1e",
"D1_Htr4",
"D1_Htr4",
"D1_Htr4",
"D1_Htr4",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D1_Ebf1",
"D2_Scube1",
"D2_Scube1",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a",
"D2_Stk32a"
)

annotation <- data.frame(cluster_group)
rownames(annotation) <- colnames(res)
annotation$cluster_group <- as.factor(annotation$cluster_group)
levels(annotation$cluster_group)
ann_colors<-list(cluster_group=c("D1_Ppm1e"="#4db81c", "D1_Htr4"="#a823d9", "D1_Ebf1"="#2399d9", "D2_Scube1"="#cf1b1b", "D2_Stk32a"="#ebbb1e"))

p7 <- as.ggplot(pheatmap(res,
			  cellheight=10, cellwidth = 10,
			  gaps_row = c(4,8,20,22),
			  gaps_col = c(4,8,20,22),
			  annotation_col=annotation,
			  annotation_legend=FALSE,
			  annotation_names_col = FALSE,
			  annotation_row=annotation,
			  annotation_names_row = FALSE,
			  annotation_colors = ann_colors,
			  fontsize = 8,
			  treeheight_row = 0, treeheight_col = 0,
			  cluster_rows = FALSE,
			  cluster_cols = FALSE,
			  angle_col=90
			  ))

layout <- "
AABBCCHHHH
AABBCCHHHH
EEFFGGHHHH
EEFFGGHHHH
"

p <- p1 + p3 + p4 + p2 + p5 + p6 + p7

p <- p +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')


ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_replication_mapquery_fig.png", width = 20, height = 7, dpi=300)
ggsave(p, filename="/home/crist/nac_rn7/nac_rn7_replication_mapquery_fig.eps", width = 20, height = 7)

sessionInfo()
