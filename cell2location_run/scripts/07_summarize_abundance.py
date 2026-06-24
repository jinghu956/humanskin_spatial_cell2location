import pandas as pd
import scanpy as sc
import os

result_h5ad = "cell2location_run/results/spatial_cell2location_broad_adult.h5ad"
abundance_csv = "cell2location_run/results/cell_abundance_broad_adult.csv"

out_clean = "cell2location_run/results/cell_abundance_broad_adult_clean.csv"
out_with_coords = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
out_summary = "cell2location_run/results/cell_abundance_broad_adult_summary.csv"
out_dominant = "cell2location_run/results/spot_dominant_celltype.csv"

adata = sc.read_h5ad(result_h5ad)
abund = pd.read_csv(abundance_csv, index_col=0)

print("Original abundance shape:", abund.shape)
print("Original columns:")
print(abund.columns.tolist())

# 清理列名
clean_cols = []
for c in abund.columns:
    c2 = c.replace("q05cell_abundance_w_sf_", "")
    c2 = c2.replace("q05_cell_abundance_w_sf_", "")
    clean_cols.append(c2)

abund.columns = clean_cols

print("\nClean columns:")
print(abund.columns.tolist())

# 保存简洁版 abundance
abund.to_csv(out_clean)

# 合并坐标
coords = adata.obs[["spot_id", "barcode", "x", "y"]].copy()
abund_with_coords = coords.join(abund)
abund_with_coords.to_csv(out_with_coords)

# 每类细胞全局统计
summary = pd.DataFrame({
    "mean_abundance": abund.mean(axis=0),
    "median_abundance": abund.median(axis=0),
    "sum_abundance": abund.sum(axis=0),
    "max_abundance": abund.max(axis=0),
    "nonzero_spots": (abund > 0).sum(axis=0),
})
summary = summary.sort_values("sum_abundance", ascending=False)
summary.to_csv(out_summary)

# 每个 spot 的主导细胞类型
dominant = pd.DataFrame(index=abund.index)
dominant["dominant_celltype"] = abund.idxmax(axis=1)
dominant["dominant_abundance"] = abund.max(axis=1)
dominant = coords.join(dominant)
dominant.to_csv(out_dominant)

print("\nSaved:")
print(out_clean)
print(out_with_coords)
print(out_summary)
print(out_dominant)

print("\nCell type abundance summary:")
print(summary)

print("\nDominant cell type counts:")
print(dominant["dominant_celltype"].value_counts())
