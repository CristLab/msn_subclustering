source("/home/crist/OUD_10x/const_plot.R")


library(future)
library(dplyr)
library(magrittr)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(scrattch.hicat)


#Modified from constellation plot code from Allen Brain and Carmen Sandoval (https://github.com/carmensandoval/singlecell-neocortex-arealization)

oudfull <- readRDS("/home/crist/nac_rn7/nac_rn7_subcluster.rds")
oudfull


head(oudfull@meta.data)




umap.coord <- as.data.frame(oudfull[["umap"]]@cell.embeddings)

cells.cl.df <- oudfull@meta.data

head(cells.cl.df)


cell_names <- rownames(cells.cl.df)
cl <- cells.cl.df[["subclusters"]] %>% as.factor %>% magrittr::set_names(value = cell_names)

#Original use of get_cl_df from scrattch.hicat fails to add the size to the output dataframe
#Not clear why. Possible version differences in some packages or formatting issues
#This modified version of the get_cl_df function code works
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", 
                                   "#7F0000"))
  
  cl.df <- data.frame(cluster_label = sort(unique(cl)))
  cl.size <- table(cl)
  cl.df$cluster_id <- 1:nrow(cl.df)
  cl.df$cluster_color <-jet.colors(nrow(cl.df))
  cl.df$size <- as.numeric(cl.size[cl.df$cluster_label])
  row.names(cl.df) <- cl.df$cluster_label


# cl.df <- get_cl_df(cl)

typeof(cl.df)
head(cl.df)


#They originally used an underscore between clade and area in the cluster_label format
#Could reformat cluster labels in original file to reuse their code for clade
cl.df$clade <- cl.df$cluster_label

cells.cl.df <- cells.cl.df %>% dplyr::rename(cluster_label = subclusters) %>%
                        left_join(cl.df, by = "cluster_label") %>%
						mutate(cluster_label = as.factor(cluster_label))

rd.cl.center <- get_RD_cl_center(rd.dat = umap.coord, cl)

rd.cl.center %<>% 
  as.data.frame %>% 
  set_names(c("x", "y")) %>%
  tibble::add_column(cl = cl.df$cluster_id, .before = "x") %>%
  tibble::rownames_to_column("cluster_label")

cl.center.df <- left_join(rd.cl.center, cl.df,
                          by = c("cluster_label")) 


cl.center.df <- cl.center.df %>% dplyr::rename(cluster_size = size)

head(cl.center.df)


knn.cl <- get_knn_graph(rd.dat = umap.coord, 
                        cl.df =  cl.df,
						cl = cl,						
                        knn.outlier.th = 2, 
                        outlier.frac.t = 0.5)

#Sets connection cutoffs for printing edges on final plot
filterKNN <- function(knn.cl, frac.th = 0.02) {
knn.cl.df.filter <- knn.cl$knn.cl.df %>% dplyr::filter(frac >= frac.th) %>% 
  mutate(cl.from = as.numeric(cl.from), cl.to = as.numeric(cl.to))
}

knn.cl.df.filter <- filterKNN(knn.cl = knn.cl)

typeof(knn.cl.df.filter)

#Without this the size and cluster_size columns are both empty
# cl.center.df[["cluster_size"]] <- c(
# 1491,
# 989,
# 909,
# 868,
# 891,
# 953,
# 360,
# 37,
# 529,
# 1794,
# 1068,
# 2075,
# 2326,
# 2571,
# 2515,
# 2693,
# 2305,
# 1873,
# 1033,
# 1495,
# 996,
# 1203,
# 1138,
# 948,
# 416,
# 1025,
# 547,
# 1677,
# 1025,
# 910,
# 713,
# 980,
# 837,
# 1681,
# 1336,
# 1925,
# 1459)

cl.center.df[["cluster_color"]] <- c(
'#154360',
'#1F618D',
'#2471A3',
'#2980B9',
'#5499C7',
'#21618C',
'#2874A6',
'#2E86C1',
'#3498DB',
'#5DADE2',
'#85C1E9',
'#7FB3D5',
'#884EA0',
'#A569BD',
'#D7BDE2',
'#6C3483',
'#2ECC71',
'#1E8449',
'#1D8348',
'#0E6655',
'#cf1b1b',
'#f54040',
'#F5B041',
'#DC7633',
'#922B21',
'#AF601A',
'#CA6F1E',
'#E67E22',
'#C0392B',
'#EB984E',
'#D68910',
'#B9770E',
'#873600',
'#A04000'
)

cl.center.df[["clade_id"]] <- c(
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Ebf1',
'D1_Htr4',
'D1_Htr4',
'D1_Htr4',
'D1_Htr4',
'D1_Ppm1e',
'D1_Ppm1e',
'D1_Ppm1e',
'D1_Ppm1e',
'D2_Scube1',
'D2_Scube1',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a',
'D2_Stk32a'
)

cl.center.df[["clade_color"]] <- c(
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#D6EAF8',
'#E8DAEF',
'#E8DAEF',
'#E8DAEF',
'#E8DAEF',
'#D5F5E3',
'#D5F5E3',
'#D5F5E3',
'#D5F5E3',
'#F7DC6F',
'#F7DC6F',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857',
'#F6B857'
)

plot.hull <- c(
'D1_Ebf1',
'D1_Htr4',
'D1_Ppm1e',
'D2_Scube1',
'D2_Stk32a'
)

plot_constellation_fill(knn.cl.df = knn.cl.df.filter, 
                              cl.center.df = cl.center.df, 
                              out.dir = "/home/crist/nac_rn7/",
                              node.label = "cluster_label",
                              exxageration = 0.5, curved = TRUE, 
							  plot.hull = plot.hull, 
                              # plot.height = 15, plot.width = 8,
                              node.dodge = TRUE,
                              label.size = 3, max_size = 10,
							  label_repel=TRUE)
							  


