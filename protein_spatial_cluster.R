# Ensure the Python environment for Giotto has been installed
genv_exists <- Giotto::checkGiottoEnvironment()

if(!"pak" %in% installed.packages()) {
  install.packages("pak")
}

pak::pkg_install("giotto-suite/Giotto")

# python environment install
library(Giotto)
#installGiottoEnvironment()

# read the data----------
data_path <- "/Users/xixigong/PycharmProjects/LIANA/human_skin/input/" #输入目录

x <- data.table::fread(paste0(data_path, "GSM6578065_humanskin_RNA.tsv"))

spatial_coords <- data.frame(cell_ID = x$X)
spatial_coords <- cbind(spatial_coords,
                        stringr::str_split_fixed(spatial_coords$cell_ID, 
                                                 pattern = "x",
                                                 n = 2))
colnames(spatial_coords)[2:3] = c("sdimx", "sdimy")
spatial_coords$sdimx <- as.integer(spatial_coords$sdimx)
spatial_coords$sdimy <- as.integer(spatial_coords$sdimy)
spatial_coords$sdimy <- spatial_coords$sdimy*(-1)

# RNA and protein
rna_matrix <- data.table::fread(paste0(data_path, "GSM6578065_humanskin_RNA.tsv"))

rna_matrix <- rna_matrix[rna_matrix$X %in% spatial_coords$cell_ID,]

rna_matrix <- rna_matrix[match(spatial_coords$cell_ID, rna_matrix$X),]

rna_matrix <- t(rna_matrix[,-1])

colnames(rna_matrix) <- spatial_coords$cell_ID

protein_matrix <- data.table::fread(paste0(data_path, "GSM6578074_humanskin_protein.tsv"))

protein_matrix <- protein_matrix[protein_matrix$X %in% spatial_coords$cell_ID,]

protein_matrix <- protein_matrix[match(spatial_coords$cell_ID, protein_matrix$X),]

protein_matrix <- t(protein_matrix[,-1])

colnames(protein_matrix) <- spatial_coords$cell_ID

# create giotto object---------

library(Giotto)

results_folder <- "~/PycharmProjects/giottosuite" # save picture to this directory

instructions <- createGiottoInstructions(save_plot = TRUE,
                                         save_dir = results_folder,
                                         show_plot = FALSE,
                                         return_plot = FALSE)


my_giotto_object <- createGiottoObject(expression = list(rna = list(raw = rna_matrix),
                                                         protein = list(raw = protein_matrix)),
                                       spatial_locs = spatial_coords,
                                       instructions = instructions)


#
img_path <- '/Users/xixigong/PycharmProjects/LIANA/Tissue_img/human_jpg/skin.jpg'

#install.packages("magick")
library(magick)
mg_img <- image_read(img_path)
my_giotto_image <- createGiottoImage(gobject = my_giotto_object,
                                     mg_object = mg_img, do_manual_adj = TRUE,
                                     scale_factor  = 0.5,negative_y    = TRUE)
my_giotto_object <- addGiottoImage(gobject = my_giotto_object, images = list(my_giotto_image), spat_loc_name = "raw")


#
spatPlot2D(my_giotto_object,
           point_size = 3.5)



# Precessing--------
# filtering
# RNA
my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                                 spat_unit = "cell",
                                 feat_type = "rna",
                                 expression_threshold = 1,
                                 feat_det_in_min_cells = 1,
                                 min_det_feats_per_cell = 1)

# protein
# Protein
my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                                 spat_unit = "cell",
                                 feat_type = "protein",
                                 expression_threshold = 1,
                                 feat_det_in_min_cells = 1,
                                 min_det_feats_per_cell = 1)

#normalization
# RNA
# RNA
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    norm_methods = "standard",
                                    scalefactor = 10000,
                                    verbose = TRUE)
# Protein
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "protein",
                                    scalefactor = 6000,
                                    verbose = TRUE)

#add cell and feature statistics(optional)
# RNA
my_giotto_object <- addStatistics(gobject = my_giotto_object,
                                  spat_unit = "cell",
                                  feat_type = "rna")
# Protein
my_giotto_object <- addStatistics(gobject = my_giotto_object,
                                  spat_unit = "cell",
                                  feat_type = "protein",
                                  expression_values = "normalized")

# Dimension Reduction--------
#PCA
# RNA
my_giotto_object <- runPCA(gobject = my_giotto_object,
                           spat_unit = "cell",
                           feat_type = "rna",
                           expression_values = "normalized",
                           reduction = "cells",
                           name = "rna.pca")
#install.packages(c("FactoMineR"))
# Protein
my_giotto_object <- runPCA(gobject = my_giotto_object,
                           spat_unit = "cell",
                           feat_type = "protein",
                           expression_values = "normalized",
                           scale_unit = TRUE,
                           center = FALSE,
                           method = "factominer")

# clustering -------
# RNA
my_giotto_object <- runUMAP(gobject = my_giotto_object,
                            spat_unit = "cell",
                            feat_type = "rna",
                            expression_values = "normalized",
                            reduction = "cells",
                            dimensions_to_use = 1:10,
                            dim_reduction_name = "rna.pca")
# Protein
my_giotto_object <- runUMAP(gobject = my_giotto_object,
                            spat_unit = "cell",
                            feat_type = "protein",
                            expression_values = "normalized",
                            dimensions_to_use = 1:10)


# create nearest network
# RNA
my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "rna",
                                         type = "sNN",
                                         dim_reduction_to_use = "pca",
                                         dim_reduction_name = "rna.pca",
                                         dimensions_to_use = 1:10,
                                         k = 20)
# Protein
my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "protein",
                                         type = "sNN",
                                         name = "protein_sNN.pca",
                                         dimensions_to_use = 1:10,
                                         k = 20)

# find leiden cluster
# RNA
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    nn_network_to_use = "sNN",
                                    name = "leiden_clus",
                                    resolution = 1)
# Protein
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "protein",
                                    nn_network_to_use = "sNN",
                                    network_name = "protein_sNN.pca",
                                    name = "leiden_clus",
                                    resolution = 1)


# Plot PCA
# RNA
plotPCA(gobject = my_giotto_object,
        spat_unit = "cell",
        feat_type = "rna",
        dim_reduction_name = "rna.pca",
        cell_color = "leiden_clus",
        title = "RNA PCA")
# Protein
plotPCA(gobject = my_giotto_object,
        spat_unit = "cell",
        feat_type = "protein",
        dim_reduction_name = "protein.pca",
        cell_color = "leiden_clus",
        title = "Protein PCA")

# PlotUMAP
# RNA
plotUMAP(gobject = my_giotto_object,
         spat_unit = "cell",
         feat_type = "rna",
         cell_color = "leiden_clus",
         point_size = 2,
         title = "RNA Uniform Manifold Approximation & Projection (UMAP)",
         axis_title = 12,
         axis_text = 10 )
# Protein
plotUMAP(gobject = my_giotto_object,
         spat_unit = "cell",
         feat_type = "protein",
         cell_color = "leiden_clus",
         dim_reduction_name = "protein.umap",
         point_size = 2,
         title = "Protein Uniform Manifold Approximation & Projection (UMAP)",
         axis_title = 12,
         axis_text = 10 )
#plot spatial locations by cluster
# RNA
spatPlot2D(my_giotto_object,
           show_image = TRUE,
           point_size = 3.5,
           cell_color = "leiden_clus",
           title = "RNA Leiden clustering")
# Protein
spatPlot2D(my_giotto_object,
           spat_unit = "cell",
           feat_type = "protein",
           cell_color = "leiden_clus",
           point_size = 3.5,
           show_image = TRUE,
           title = "Protein Leiden clustering")

# multi-omics integration--------
# knn to create nearest network
my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "rna",
                                         type = "kNN",
                                         dim_reduction_name = "rna.pca",
                                         name = "rna_kNN.pca",
                                         dimensions_to_use = 1:10,
                                         k = 20)
my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "protein",
                                         type = "kNN",
                                         name = "protein_kNN.pca",
                                         dimensions_to_use = 1:10,
                                         k = 20)
# run WNN
my_giotto_object <- runWNN(
  gobject   = my_giotto_object,
  feat_type = c("rna", "protein"),
  k         = 20
)

# integrated umap
my_giotto_object <- runIntegratedUMAP(my_giotto_object,
                                      feat_types = c("rna", "protein"))

# integrated leiden clusters
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    nn_network_to_use = "kNN",
                                    network_name = "integrated_kNN",
                                    name = "integrated_leiden_clus",
                                    resolution = 0.5)

plotUMAP(gobject = my_giotto_object,
         spat_unit = "cell",
         feat_type = "rna",
         cell_color = "integrated_leiden_clus",
         dim_reduction_name = "integrated.umap",
         point_size = 2.5,
         title = "Integrated UMAP using Integrated Leiden clusters",
         axis_title = 12,
         axis_text = 10)
# integrated spatial locations
spatPlot2D(my_giotto_object,
           spat_unit = "cell",
           feat_type = "rna",
           cell_color = "integrated_leiden_clus",
           point_size = 3.5,
           show_image = TRUE,
           title = "Integrated Leiden clustering")

# spatially variable genes-----------
my_giotto_object <- createSpatialNetwork(gobject = my_giotto_object,
                                         method = "kNN", 
                                         k = 6,
                                         maximum_distance_knn = 5, #distance cuttof for nearest neighbors to consider for kNN network
                                         name = "spatial_network")

rank_spatialfeats <- binSpect(my_giotto_object, 
                              bin_method = "rank",
                              calc_hub = TRUE, 
                              hub_min_int = 5,
                              spatial_network_name = "spatial_network")

spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = rank_spatialfeats$feats[1:6], 
               cow_n_col = 2, 
               point_size = 1.5)

# spatially corelated genes
# 3.1 cluster the top 500 spatial genes into 20 clusters
ext_spatial_genes <- rank_spatialfeats[1:500,]$feats

# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT <- detectSpatialCorFeats(my_giotto_object,
                                          method = "network",
                                          spatial_network_name = "spatial_network",
                                          subset_feats = ext_spatial_genes)


# 3.3 identify potential spatial co-expression
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT, 
                                           name = "spat_netw_clus", 
                                           k = 3)
# install.packages(c("ComplexHeatmap"))
# visualize clusters: failed because ComplexHeatmap和R4.5不兼容
heatmSpatialCorFeats(my_giotto_object,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL),
                     save_param = list(base_height = 6, 
                                       base_width = 8, 
                                       units = "cm"))

# Metagenes/co-expression modules----------
# 3.4 create metagenes / co-expression modules
cluster_genes <- getBalancedSpatCoexpressionFeats(spat_cor_netw_DT, 
                                                  maximum = 30)

my_giotto_object <- createMetafeats(my_giotto_object, 
                                    feat_clusters = cluster_genes, 
                                    name = "cluster_metagene")

spatCellPlot(my_giotto_object,
             spat_enr_names = "cluster_metagene",
             cell_annotation_values = as.character(c(1:7)),
             point_size = 1, 
             cow_n_col = 3)


# Spatially informed clustering--------
my_spatial_genes = names(cluster_genes)

my_giotto_object <- runPCA(gobject = my_giotto_object,
                           feats_to_use = my_spatial_genes,
                           name = "custom_pca")

my_giotto_object <- runUMAP(my_giotto_object,
                            dim_reduction_name = "custom_pca",
                            dimensions_to_use = 1:20,
                            name = "custom_umap")

my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         dim_reduction_name = "custom_pca",
                                         dimensions_to_use = 1:20, 
                                         k = 3,
                                         name = "custom_NN")

my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    network_name = "custom_NN",
                                    resolution = 0.05, 
                                    n_iterations = 1000,
                                    name = "custom_leiden")
spatPlot2D(my_giotto_object,
           show_image = FALSE,
           cell_color = "custom_leiden",
           point_size = 3.5,
           background_color = "black",
           title = "Spatially informed clustering")

spatPlot2D(my_giotto_object,
           show_image = FALSE,
           cell_color = "custom_leiden", 
           cell_color_code = c(
             "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
             "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
             "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5"
           ),
           point_size = 3.5,
           background_color = "black",
           title = "Spatially informed clustering")

spatPlot2D(my_giotto_object,
           show_image = FALSE,
           cell_color = "custom_leiden", 
           # 👇 正好 15 个颜色，适配你的15个cluster
           cell_color_code = c(
             "#eb4034",    # 1
             "#5877e8",    # 2
             "#ebd834",    # 3
             "#9beb34",    # 4
             "#6fab6a",    # 5
             "#24703f",    # 6
             "#58e8cb",    # 7
             "#58d0e8",    # 8
             "#eb8f34",    # 9
             "#7f58e8",    #10
             "#d758e8",    #11
             "#e85892",    #12
             "#FF69B4",    #13 新增
             "#ADFF2F",    #14 新增
             "#FF4500"     #15 新增
           ),
           point_size = 3.5,
           background_color = "black",
           title = "Spatially informed clustering")


# Spatially variable proteins----------
# 用 protein 层创建空间网络（如果还没创建过，复用之前的也可以）
# create nearest network

# binSpect 指定 feat_type = "protein"
rank_spatialfeats_protein <- binSpect(my_giotto_object,
                                      spat_unit = "cell",
                                      feat_type = "protein",   # ← 关键改动
                                      bin_method = "rank", #method to binarize gene expression
                                      calc_hub = TRUE,
                                      hub_min_int = 5,
                                      spatial_network_name = "spatial_network")

# Spatially correlated proteins----------
# 检查 protein 层是否有 normalized 表达矩阵
my_giotto_object@expression$cell$protein
# 应该有 normalized 或 scaled 层

# protein panel 一般只有几十到几百个，不需要截取500个
# 根据你实际的 protein 数量调整，比如全用或取 top N
ext_spatial_proteins <- rank_spatialfeats_protein$feats  # 全用，或 [1:N]

spat_cor_netw_DT_protein <- detectSpatialCorFeats( #pairwise distances between genes
  my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",            # ← 关键改动
  method = "network",
  spatial_network_name = "spatial_network",
  subset_feats = ext_spatial_proteins
)

spat_cor_netw_DT_protein <- clusterSpatialCorFeats(
  spat_cor_netw_DT_protein,
  name = "spat_netw_clus_protein",  # ← 换个名字，避免覆盖 RNA 结果
  k = 3                             # protein 类型少，k 可以适当调小
)

# Metagenes（Metaproteins）----------
cluster_proteins <- getBalancedSpatCoexpressionFeats(
  spat_cor_netw_DT_protein,
  maximum = 30
)

my_giotto_object <- createMetafeats(
  my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",               # ← 关键改动
  feat_clusters = cluster_proteins,
  name = "cluster_metaprotein"         # ← 换名字
)
# 看看 enrichment 结果存在哪里
slot(my_giotto_object, "spatial_enrichment")
# 确认 metafeats 存在哪里
showGiottoSpatEnrichments(my_giotto_object)

spatCellPlot(my_giotto_object,
             feat_type = 'protein',# must set,feature type (e.g. "rna", "dna", "protein")
             spat_enr_names = "cluster_metaprotein",
             cell_annotation_values = as.character(c(1:3)),  # 对应 k=3
             point_size = 1,
             cow_n_col = 3)


# Spatially informed clustering (protein)----------
my_spatial_proteins <- names(cluster_proteins)

my_giotto_object <- runPCA(gobject = my_giotto_object,
                           spat_unit = "cell",
                           feat_type = "protein",        # ← 关键改动
                           feats_to_use = my_spatial_proteins,
                           name = "custom_pca_protein")  # ← 换名字

my_giotto_object <- runUMAP(my_giotto_object,
                            spat_unit = "cell",
                            feat_type = "protein",       # ← 关键改动
                            dim_reduction_name = "custom_pca_protein",
                            dimensions_to_use = 1:10,    # protein 维度通常比 RNA 少
                            name = "custom_umap_protein")

my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "protein",      # ← 关键改动
                                         type = 'kNN', # 明确指定
                                         dim_reduction_name = "custom_pca_protein", # 用proteinPCA
                                         dimensions_to_use = 1:10,
                                         k = 3,
                                         name = "custom_NN_protein_kNN") # ← 换名字

my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "protein",           # ← 关键改动
                                    nn_network_to_use = "kNN",
                                    network_name = "custom_NN_protein_kNN",
                                    resolution = 0.5,
                                    n_iterations = 1000,
                                    name = "custom_leiden_protein")  # ← 换名字

spatPlot2D(my_giotto_object,
           show_image = TRUE,
           spat_unit = "cell",
           feat_type = "protein",
           cell_color = "custom_leiden_protein",
           point_size = 3.5,
           background_color = "black",
           title = "Spatially informed clustering (Protein)")

spatPlot2D(my_giotto_object,
           #show_image = TRUE,
           spat_unit = "cell",
           feat_type = "protein",
           cell_color = "custom_leiden_protein",
           point_size = 3.5,
           #background_color = "black",
           title = "Spatially informed clustering (Protein)")


# 对比两个结果
# 并排对比
library(ggplot2)
p1 <- spatPlot2D(my_giotto_object, cell_color = "custom_leiden",
                 title = "RNA-based", return_plot = TRUE)
p2 <- spatPlot2D(my_giotto_object, cell_color = "custom_leiden_protein",
                 title = "Protein-based", return_plot = TRUE)
cowplot::plot_grid(p1, p2)

# save marker protein
cluster1marker = names(cluster_proteins[cluster_proteins == 1])
cluster2marker = names(cluster_proteins[cluster_proteins == 2])
cluster3marker = names(cluster_proteins[cluster_proteins == 3])

cluster1markerc <- str_sub(cluster1marker, 1, -16)
cluster2markerc <- str_sub(cluster2marker, 1, -16)
cluster3markerc <- str_sub(cluster3marker, 1, -16)

write.table(cluster3markerc,'cluster3.txt',row.names = FALSE,quote = FALSE)
write.table(cluster2markerc,'cluster2.txt',row.names = FALSE,quote = FALSE)
write.table(cluster1markerc,'cluster1.txt',row.names = FALSE,quote = FALSE)

df_cluster <- data.frame(cluster1=cluster1markerc,
                        cluster2=cluster2markerc,
                        cluster3=cluster3markerc)





# show marker gene
spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = cluster1marker[1:6] , 
               feat_type = 'protein',
               cow_n_col = 2, 
               point_size = 1.5)   
spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = cluster1marker[7:12] , 
               feat_type = 'protein',
               cow_n_col = 2, 
               point_size = 1.5)
spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = cluster1marker[13:18] , 
               feat_type = 'protein',
               cow_n_col = 2, 
               point_size = 1.5)
spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = cluster1marker[19:24] , 
               feat_type = 'protein',
               cow_n_col = 2, 
               point_size = 1.5)
spatFeatPlot2D(my_giotto_object, 
               expression_values = "scaled",
               feats = cluster1marker[25:30] , 
               feat_type = 'protein',
               cow_n_col = 2, 
               point_size = 1.5)



for (k in 1:5){
  spatFeatPlot2D(
    my_giotto_object,
    expression_values = "scaled",
    feats = cluster2marker[(6*k - 5):(6*k)],
    feat_type = "protein",
    cow_n_col = 2,
    point_size = 1.5
  )
}

for (k in 1:3){
  idx <- if (k == 3) (6*k - 5):17 else (6*k - 5):(6*k)
  
  spatFeatPlot2D(
      my_giotto_object,
      expression_values = "scaled",
      feats = cluster3marker[idx],
      feat_type = "protein",
      cow_n_col = 2,
      point_size = 1.5 )
}

spatFeatPlot2D(my_giotto_object,
               feat_type = "protein",
               expression_values = "normalized",
               feats = c("CD207..Langerin..",   # 表皮层
                         "CD31.",               # 真皮血管
                         "Podoplanin.",         # 淋巴管
                         "CD206..MMR..",        # 真皮巨噬细胞
                         "CD3.",                # T细胞
                         "CD45."),              # 泛免疫
               cow_n_col = 2,
               point_size = 1.5)




# 不同custom_leiden_protein ----------
markers_gini <- findMarkers_one_vs_all(
  gobject = my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",
  method = "gini",
  expression_values = "normalized",
  cluster_column = "custom_leiden_protein",
  min_expr_gini_score = 0.5,
  min_det_gini_score = 0.5,
  min_feats = 1
)

head(markers_gini)

table(markers_gini$cluster)

#
# 提取每个簇的 top markers
topgenes_protein <- markers_gini[, head(.SD, 3), by = "cluster"]$feats

plotMetaDataHeatmap(
  gobject = my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",
  expression_values = "normalized",
  metadata_cols = "custom_leiden_protein",
  selected_feats = unique(topgenes_protein)
)

# violin plot
violinPlot(
  gobject = my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",
  expression_values = "normalized",
  feats = topgenes_protein[1:6],          # 选几个看
  cluster_column = "custom_leiden_protein"
)


# UMAP 上展示某个protein的表达
dimFeatPlot2D(
  gobject = my_giotto_object,
  spat_unit = "cell",
  feat_type = "protein",
  expression_values = "normalized",
  feats = topgenes_protein[1:4],
  dim_reduction_name = "custom_umap_protein",
  cow_n_col = 2,
  point_size = 1
)


# 保存差异 protein 结果
write.csv(markers_gini, 
          file = "protein_markers_per_cluster.csv", 
          row.names = FALSE)

#------
spatFeatPlot2D(my_giotto_object,
               feat_type = "protein",
               expression_values = "normalized",
               feats = c("CD207..Langerin..CATTCTTCACGGGAT",   # 表皮层
                         "CD317..BST2_.Tetherin..AAGAGCCGTTGTGAA",     # 真皮血管
                         "Podoplanin.GGTTACTCGTTGTGT",         # 淋巴管
                         "CD206..MMR..TCAGAACGTCTAACT",        # 真皮巨噬细胞
                         "CD3.CTCATTGTAACTCCT",                # T细胞
                         "CD45.TCCCTTGCGATTTAC"),              # 泛免疫
               cow_n_col = 2,
               point_size = 1.5)



