library(future)
library(dplyr)
library(Seurat)
library(SoupX)
library(DoubletFinder)
library(ggplot2)


input.list <- scan("/home/crist/nac_rn7/nac_rn7_input.txt", what = list("numeric","character","numeric","character","character","character","character","character","logical","character","numeric","numeric"))

sample.size <- input.list[[1]]
species <- input.list[[2]]
mito.cutoff <- as.numeric(input.list[[3]])
path.to.files <- input.list[[4]]
path.to.save <- input.list[[5]]
metadata.path <- input.list[[6]]
project <- input.list[[7]]
dot.features <- unlist(strsplit(input.list[[8]],","))
use.reference <- input.list[[9]]
references <- as.numeric(unlist(strsplit(input.list[[10]],",")))
cluster.res1 <- as.numeric(input.list[[11]])
cluster.res2 <- as.numeric(input.list[[12]])

#Create list of sample names. IS THIS NEEDED ANYMORE?
sample.names <- list.files(path.to.files)

#Create list of directory name for each sample. Print list to ensure order is correct
sample.dirs <- list.files(path.to.files, full.names=TRUE)
sample.dirs

print("Reading matrices") 
data.list <- lapply(X   = sample.dirs,
                         FUN = function(x){
                            Read10X(data.dir=paste0(x,"/outs/filtered_feature_bc_matrix"))
                         })

print("Making Seurat objects") 
seurat.list <- lapply(X   = data.list,
                         FUN = function(x){
                            CreateSeuratObject(counts = x, min.cells=1, min.features=1)
                         })


for (i in seq(1,sample.size)){
	print(seurat.list[[i]])
}

#Calculate the mitochondrial fractions
print("Calculating mitochondrial fraction") 
if (species == "human") {
seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^MT-", 
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
} else if (species == "mouse") {
seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^mt-", 
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
} else if (species == "rat") {
seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^Mt-",
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^AY[0-9]",
                                                col.name = "percent.ay",
                                                assay    = "RNA")
                         })

lapply(X   = seurat.list,
                         FUN = function(x){
						   x@meta.data$percent.mt <- as.numeric(x@meta.data$percent.mt) + as.numeric(x@meta.data$percent.ay)
                        })
} else {
print("Species not found")
stop()
}

#Basic plots						 
lapply(X   = seurat.list,
			FUN = function(x){
			print(head(x@meta.data$percent.mt))
            VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
            })




#Basic QC					 
seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                           x <- subset(x, subset = nFeature_RNA > 200 & percent.mt < mito.cutoff)
                         })

for (i in seq(1,sample.size)){
	print(seurat.list[[i]])
}

#Post QC plots						 
lapply(X   = seurat.list,
			FUN = function(x){
            VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
            })


#Read in metadata file
metadata <- read.csv(file = metadata.path ,header=TRUE)

#Add metadata characteristics
for (column in colnames(metadata)){
	for (i in seq(1,sample.size)){
							seurat.list[[i]] <- AddMetaData(seurat.list[[i]], metadata[[as.character(column)]][[i]], col.name = as.character(column))
                         }
}


#Use SCTransform to integrate
plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2, future.seed=TRUE)

seurat.list <- lapply(X   = seurat.list,
                         FUN = function(x){
                            SCTransform(x, vst.flavor = "v2", vars.to.regress = c("nCount_RNA"), verbose = FALSE) %>%
							RunPCA(npcs = 20, verbose = FALSE)
                         })

features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)


if (use.reference) {
anchors <- FindIntegrationAnchors(object.list = seurat.list, 
									reference = references, 
									normalization.method = "SCT",
     								anchor.features = features)
} else {
anchors <- FindIntegrationAnchors(object.list = seurat.list,  
									normalization.method = "SCT",
     								anchor.features = features)
}

sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")


sct <- RunPCA(sct, npcs = 20, verbose = FALSE)
sct <- RunUMAP(sct, reduction = "pca", dims = 1:20, verbose = FALSE)
sct <- FindNeighbors(sct, reduction = "pca", dims = 1:20)
sct <- FindClusters(sct, resolution = cluster.res1)
sct

#Print tables for clusters by metadata columns
for (column in colnames(metadata)){
	print(table(Idents(sct),sct@meta.data[[as.character(column)]]))
}

#Plot starting UMAP
DimPlot(sct, reduction = "umap", group.by='ident') + NoLegend()

#Dot plot for major markers using RNA assay
DefaultAssay(sct) <- "RNA"

p <- DotPlot(sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename=paste0(path.to.save,"dotplot_",project,"sct_1.png"))

saveRDS(sct, file = paste0(path.to.save,project,"_sct_working.rds"))

#Cell metadata column
sct$cell <- row.names(sct@meta.data)

#Strip underscore and number off to get original cell ID so it matches raw/filtered matrices
sct$orig_cellid <- as.character(lapply(strsplit(x = sct$cell,split = "_"),"[",1))

#Create an empty vector that will hold the matrices and single GEM objects
adjusted_matrices <- vector(mode = "list",length = sample.size)
seurat_objects    <- vector(mode = "list",length = sample.size)


#SoupX pipeline
for (i in seq(1,sample.size)){
	filtered_mat <- Seurat::Read10X(file.path(paste0(sample.dirs[i],"/outs/filtered_feature_bc_matrix/")))
	raw_mat <- Seurat::Read10X(file.path(paste0(sample.dirs[i],"/outs/raw_feature_bc_matrix/")))
	sc  <- SoupChannel(tod = raw_mat, toc = filtered_mat)
	#Add Metadata info that includes cell and celltype
	sc$metaData$cell <- row.names(sc$metaData)
	sc$metaData$celltype <- apply(X = sc$metaData,MARGIN = 1,function(x){
						as.character(subset(sct@meta.data, subset=(orig_cellid == x["cell"] & sample_num == i))$seurat_clusters)
						})
	sc <- setClusters(sc,clusters = as.character(sc$metaData$celltype))
	sc  <- autoEstCont(sc,forceAccept = TRUE)
	out <- adjustCounts(sc,roundToInt = TRUE)
	adjusted_matrices[[i]] <- out
	meta <-  sct@meta.data[which(sct$orig_cellid %in% colnames(adjusted_matrices[[i]]) & sct$sample_num == i),]
	row.names(meta) <- meta$orig_cellid
	seurat_objects[[i]] <- CreateSeuratObject(adjusted_matrices[[i]],
                                                       meta.data = meta)
													   
  #Remove the rows with no metadata												
  seurat_objects[[i]] <- seurat_objects[[i]][,colnames(seurat_objects[[i]]) %in% row.names(meta)]												   
														   
													   
  rm(sc,out,meta,filtered_mat,raw_mat)
  
}


soupx <- seurat_objects

#Add new mito percentage
if (species == "human") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^MT-", 
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
} else if (species == "mouse") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^mt-", 
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
} else if (species == "rat") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^Mt-",
                                                col.name = "percent.mt",
                                                assay    = "RNA")
                         })
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                pattern = "^AY[0-9]",
                                                col.name = "percent.ay",
                                                assay    = "RNA")
                         })

lapply(X   = soupx,
                         FUN = function(x){
						   x@meta.data$percent.mt <- as.numeric(x@meta.data$percent.mt) + as.numeric(x@meta.data$percent.ay)
                        })
} else {
print("Species not found")
stop()
}

#Remove mt transcripts

if (species == "human") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                            all.genes <- as.data.frame(rownames(x))
							colnames(all.genes) <- c("genes")
							nomtgenes <- all.genes[!grepl("^MT-",all.genes$genes),]
							x <- subset(x, features = nomtgenes)
                         })
} else if (species == "mouse") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                            all.genes <- as.data.frame(rownames(x))
							colnames(all.genes) <- c("genes")
							nomtgenes <- all.genes[!grepl("^mt-",all.genes$genes),]
							x <- subset(x, features = nomtgenes)
                         })
} else if (species == "rat") {
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                            all.genes <- as.data.frame(rownames(x))
							colnames(all.genes) <- c("genes")
							nomtgenes <- all.genes[!grepl("^Mt-",all.genes$genes),]
							x <- subset(x, features = nomtgenes)
                         })
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                            all.genes <- as.data.frame(rownames(x))
							colnames(all.genes) <- c("genes")
							noaygenes <- all.genes[!grepl("^AY[0-9]",all.genes$genes),]
							x <- subset(x, features = noaygenes)
                         })
} else {
print("Species not found")
stop()
}

#Plot the number of genes, UMIs, and percent.mt for every cell
lapply(X   = soupx,
       FUN = function(x){
         VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
       })


#Now subset for cells with > 200 genes, less than 5k genes, and less than 5% of reads mapping to the mitochondrial genome
soupx <- lapply(X   = soupx,
                         FUN = function(x){
                           subset(x = x, subset  =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < mito.cutoff)
                         })

#SCTransform to normalize/scale/integrate
plan("multicore",workers=16)
options(future.globals.maxSize = 1000000 * 1024^2, future.seed=TRUE)

soupx <- lapply(X   = soupx,
                         FUN = function(x){
                            SCTransform(x, vst.flavor = "v2", vars.to.regress = c("nCount_RNA"), verbose = FALSE) %>%
							RunPCA(npcs = 20, verbose = FALSE)
                         })


features <- SelectIntegrationFeatures(object.list = soupx, nfeatures = 3000)
soupx <- PrepSCTIntegration(object.list = soupx, anchor.features = features)

if (use.reference) {
soupx.anchors <- FindIntegrationAnchors(object.list = soupx, 
									reference = references, 
									normalization.method = "SCT",
     								anchor.features = features)
} else {
soupx.anchors <- FindIntegrationAnchors(object.list = soupx,  
									normalization.method = "SCT",
     								anchor.features = features)
}

soupx.sct <- IntegrateData(anchorset = soupx.anchors, normalization.method = "SCT")


soupx.sct <- RunPCA(soupx.sct, npcs = 20, verbose = FALSE)
soupx.sct <- RunUMAP(soupx.sct, reduction = "pca", dims = 1:20, verbose = FALSE)
soupx.sct <- FindNeighbors(soupx.sct, reduction = "pca", dims = 1:20)
soupx.sct <- FindClusters(soupx.sct, resolution = cluster.res2)
soupx.sct

#Violin plot of UMI for each cluster
VlnPlot(soupx.sct, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

#Print number of nuclei per sample per cluster
for (column in colnames(metadata)){
	print(table(Idents(soupx.sct),soupx.sct@meta.data[[as.character(column)]]))
}

#Plot starting UMAP
DimPlot(soupx.sct, reduction = "umap", group.by='ident') + NoLegend()

#Plot UMI per cluster
VlnPlot(soupx.sct, features = c("nCount_RNA"), pt.size=0) + theme(legend.position = 'none')

#Dot plot for major markers using RNA assay
DefaultAssay(soupx.sct) <- "RNA"

p <- DotPlot(soupx.sct, features=dot.features, cols=c("lightgrey","midnightblue"), col.min=(-1), scale.min=0) +
	guides(color = guide_colorbar(title = 'Avg. Exp.'),size = guide_legend(title = 'Pct. Exp.')) +
	theme(
		axis.text.y=element_text(hjust=0,face="bold"),
		axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold.italic"),
		axis.title.x=element_blank(),
		axis.title.y=element_blank()
	)

ggsave(p, filename=paste0(path.to.save,"dotplot_",project,"sct_soupx_1.png"))

saveRDS(soupx.sct, file = paste0(path.to.save,project,"_sct_soupx.rds"))

sessionInfo()


	