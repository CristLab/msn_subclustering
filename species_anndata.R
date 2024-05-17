library(future)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratDisk)

library(reticulate)


##Make human covariate and h5ad file

load(file='/home/crist/OUD_10x/Human_NAc/SCE_NAc-n8_tran-etal.rda')

human_nac <- as.Seurat(sce.nac.tran)
human_nac

Idents(human_nac) <- 'cellType'

saveRDS(human_nac, file = "/home/crist/nac_rn7/human_nac.rds")

head(human_nac@meta.data)

table(human_nac@meta.data$cellType)

human_msn <- subset(x = human_nac, idents = c("MSN.D1_A","MSN.D1_B","MSN.D1_C","MSN.D1_D","MSN.D1_E","MSN.D1_F","MSN.D2_A","MSN.D2_B","MSN.D2_C","MSN.D2_D"))
human_msn

human_msn@meta.data$cell <- rownames(human_msn@meta.data)

covariate <- human_msn@meta.data[,c("cell","nFeature_originalexp")]

write.table(covariate, "/home/crist/nac_rn7/scDRS/human_msn.cov", sep = "\t", row.names = FALSE)

human_msn <- DietSeurat(human_msn)
human_msn

SaveH5Seurat(human_msn, filename = "human_msn.h5Seurat")
Convert("human_msn.h5Seurat", dest = "h5ad")

#Make mouse covariate and h5ad file

mouse_msn <- readRDS("/home/crist/nac_rn7/mouse_msn_clean.rds")
mouse_msn

mouse_msn@meta.data$cell <- rownames(mouse_msn@meta.data)

head(mouse_msn@meta.data)

covariate <- mouse_msn@meta.data[,c("cell","nFeature_RNA")]

write.table(covariate, "/home/crist/nac_rn7/scDRS/mouse_msn.cov", sep = "\t", row.names = FALSE)

mouse_msn <- DietSeurat(mouse_msn)
mouse_msn

SaveH5Seurat(mouse_msn, filename = "mouse_msn.h5Seurat")
Convert("mouse_msn.h5Seurat", dest = "h5ad")

##Make rat covariate and h5ad file

rat_msn <- readRDS("/home/crist/nac_rn7/nac_rn7_subcluster.rds")
rat_msn

DefaultAssay(rat_msn) <- "RNA"

rat_msn@meta.data$cell <- rownames(rat_msn@meta.data)

covariate <- rat_msn@meta.data[,c("cell","nFeature_RNA")]

write.table(covariate, "/home/crist/nac_rn7/scDRS/rat_msn.cov", sep = "\t", row.names = FALSE)

rat_msn <- DietSeurat(rat_msn)
rat_msn

SaveH5Seurat(rat_msn, filename = "rat_msn.h5Seurat")
Convert("rat_msn.h5Seurat", dest = "h5ad")


sessionInfo()

