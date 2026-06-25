suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
})

base_dir <- "~/human_skin_for_hj"

rna_file <- file.path(base_dir, "input/GSM6578065_humanskin_RNA.tsv")

out_dir <- file.path(base_dir, "giotto_run")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. 读取 RNA 矩阵
# =========================

message("Reading RNA matrix...")

rna_dt <- fread(rna_file, data.table = FALSE)

message("RNA raw dim:")
print(dim(rna_dt))

message("RNA first columns:")
print(colnames(rna_dt)[1:8])

# 第一列是空间点位 ID，例如 15x34
spot_ids <- as.character(rna_dt[[1]])

# 删除第一列，剩下全是基因表达
rna_dt[[1]] <- NULL

rna_mat <- as.matrix(rna_dt)
rownames(rna_mat) <- spot_ids

# 当前矩阵是 spots x genes
# Giotto 需要 genes x spots
rna_mat <- t(rna_mat)

# 转成 sparse matrix
rna_mat <- Matrix(rna_mat, sparse = TRUE)

message("RNA matrix dim after transpose, genes x spots:")
print(dim(rna_mat))

# =========================
# 2. 从 spot ID 中解析空间坐标
# =========================

message("Parsing spatial coordinates from spot IDs...")

coord_split <- strsplit(colnames(rna_mat), "x")

sdimx <- as.numeric(sapply(coord_split, `[`, 1))
sdimy <- as.numeric(sapply(coord_split, `[`, 2))

if (any(is.na(sdimx)) || any(is.na(sdimy))) {
  stop("Failed to parse spatial coordinates from spot IDs.")
}

spatial_locs <- data.frame(
  cell_ID = colnames(rna_mat),
  sdimx = sdimx,
  sdimy = sdimy
)

message("Spatial locs head:")
print(head(spatial_locs))

message("Number of genes:")
print(nrow(rna_mat))

message("Number of spots:")
print(ncol(rna_mat))

# =========================
# 3. 创建 Giotto object
# =========================

message("Creating Giotto object...")

instrs <- createGiottoInstructions(
  save_dir = out_dir,
  save_plot = TRUE,
  show_plot = FALSE
)

gobj <- createGiottoObject(
  expression = rna_mat,
  spatial_locs = spatial_locs,
  instructions = instrs
)

saveRDS(gobj, file.path(out_dir, "01_giotto_object_raw.rds"))


# =========================
# 5. QC + 过滤
# =========================

message("Normalizing before statistics...")

gobj <- normalizeGiotto(gobj)

message("Adding statistics...")

gobj <- addStatistics(
  gobj,
  expression_values = "normalized"
)

# 先画 QC 图
message("Plotting QC...")


message("Available cell metadata columns after addStatistics:")

cell_meta_cols <- colnames(pDataDT(gobj))
print(cell_meta_cols)

pick_col <- function(candidates, available) {
  x <- intersect(candidates, available)
  if (length(x) == 0) {
    return(NA_character_)
  } else {
    return(x[1])
  }
}

feat_qc_col <- pick_col(
  c("nr_feats", "nr_genes", "detected_feats", "detected_genes", "n_feats", "n_genes"),
  cell_meta_cols
)

expr_qc_col <- pick_col(
  c("total_expr", "total_counts", "sum_expr", "n_counts", "total_count", "total"),
  cell_meta_cols
)

message("Feature QC column selected: ", feat_qc_col)
message("Expression QC column selected: ", expr_qc_col)

if (!is.na(feat_qc_col)) {
  png(file.path(out_dir, paste0("QC_spat_", feat_qc_col, ".png")), width = 1600, height = 1200, res = 180)
  spatPlot2D(
    gobj,
    cell_color = feat_qc_col,
    point_size = 1.5
  )
  dev.off()
} else {
  message("No feature QC column found, skipping feature QC plot.")
}

if (!is.na(expr_qc_col)) {
  png(file.path(out_dir, paste0("QC_spat_", expr_qc_col, ".png")), width = 1600, height = 1200, res = 180)
  spatPlot2D(
    gobj,
    cell_color = expr_qc_col,
    point_size = 1.5
  )
  dev.off()
} else {
  message("No expression QC column found, skipping expression QC plot.")
}

message("Filtering Giotto object...")

gobj <- filterGiotto(
  gobj,
  expression_threshold = 1,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 100
)

saveRDS(gobj, file.path(out_dir, "02_giotto_object_filtered.rds"))

# 过滤之后重新归一化
message("Re-normalizing after filtering...")

gobj <- normalizeGiotto(gobj)

message("Adding statistics after filtering...")

gobj <- addStatistics(
  gobj,
  expression_values = "normalized"
)

# =========================
# 6. 标准预处理
# =========================

message("Finding highly variable features...")

gobj <- calculateHVF(gobj)

message("Running PCA...")

gobj <- runPCA(
  gobj,
  feats_to_use = NULL,
  scale_unit = TRUE,
  center = TRUE
)

message("Running UMAP...")

gobj <- runUMAP(
  gobj,
  dimensions_to_use = 1:20
)

message("Building nearest network...")

gobj <- createNearestNetwork(
  gobj,
  dimensions_to_use = 1:20,
  k = 20
)

message("Running Leiden clustering...")

gobj <- doLeidenCluster(
  gobj,
  resolution = 0.4,
  n_iterations = 1000,
  name = "leiden_clus"
)

saveRDS(gobj, file.path(out_dir, "03_giotto_object_processed.rds"))

# =========================
# 7. 输出核心图
# =========================

message("Plotting UMAP cluster...")


png(file.path(out_dir, "UMAP_leiden_cluster.png"), width = 1600, height = 1200, res = 180)
dimPlot2D(
  gobj,
  dim_reduction_to_use = "umap",
  cell_color = "leiden_clus",
  point_size = 1.2,
  show_NN_network = FALSE
)
dev.off()

message("Plotting spatial cluster...")

png(file.path(out_dir, "Spatial_leiden_cluster.png"), width = 1600, height = 1200, res = 180)
spatPlot2D(
  gobj,
  cell_color = "leiden_clus",
  point_size = 1.5
)
dev.off()

# =========================
# 8. 画几个皮肤常用 marker
# =========================


marker_genes <- c(
  "KRT14", "KRT5", "KRT10", "KRT1",
  "COL1A1", "COL1A2", "DCN",
  "PECAM1", "VWF",
  "PTPRC", "CD3D", "MS4A1",
  "LYZ", "LST1",
  "PMEL", "MLANA"
)

current_feats <- fDataDT(gobj)$feat_ID

message("Number of features after filtering:")
print(length(current_feats))

message("Markers requested:")
print(marker_genes)

marker_genes_found <- intersect(marker_genes, current_feats)
marker_genes_missing <- setdiff(marker_genes, current_feats)

message("Marker genes found after filtering:")
print(marker_genes_found)

message("Marker genes missing after filtering:")
print(marker_genes_missing)

if (length(marker_genes_found) > 0) {
  for (gene in marker_genes_found) {
    message("Plotting marker: ", gene)

    tryCatch({
      png(
        file.path(out_dir, paste0("Spatial_", gene, ".png")),
        width = 1600,
        height = 1200,
        res = 180
      )

      spatFeatPlot2D(
        gobj,
        expression_values = "normalized",
        feats = gene,
        point_size = 1.5
      )

      dev.off()
    }, error = function(e) {
      message("Failed to plot ", gene, ": ", e$message)
      try(dev.off(), silent = TRUE)
    })
  }
} else {
  message("No marker genes found after filtering. Skipping marker plots.")
}

message("Done.")
message("Results saved to: ", out_dir)