import pandas as pd
import anndata as ad

sp_path = "input/GSM6578065_humanskin_RNA.tsv"
barcode_path = "input/spatial_barcodes_index.txt"
out_path = "cell2location_run/processed/spatial_rna.h5ad"

# RNA 矩阵：spot × gene，第一列 X 是 spot ID，例如 15x34
df = pd.read_csv(sp_path, sep="\t", index_col=0)

print("Original spatial RNA dataframe:", df.shape)

# 去掉 ERCC spike-in
keep_genes = [g for g in df.columns if not str(g).startswith("ERCC")]
df = df.loc[:, keep_genes]

print("After removing ERCC:", df.shape)

# 构建 AnnData
adata_vis = ad.AnnData(X=df.values)
adata_vis.obs_names = df.index.astype(str)
adata_vis.var_names = df.columns.astype(str)
adata_vis.var_names_make_unique()

# 读取 barcode + 坐标文件
barcodes = pd.read_csv(
    barcode_path,
    sep=r"\s+",
    header=None,
    names=["barcode", "x", "y"]
)

barcodes["spot_id"] = barcodes["x"].astype(str) + "x" + barcodes["y"].astype(str)
barcodes = barcodes.set_index("spot_id")

# 检查坐标匹配
overlap = adata_vis.obs_names.intersection(barcodes.index)
print("Spatial RNA spots:", adata_vis.n_obs)
print("Barcode file spots:", barcodes.shape[0])
print("Overlap spots:", len(overlap))

if len(overlap) != adata_vis.n_obs:
    missing = adata_vis.obs_names.difference(barcodes.index)
    print("Warning: some RNA spots are not found in barcode file.")
    print("Number missing:", len(missing))
    print("First missing examples:", missing[:10].tolist())

# 加入 barcode 和坐标
adata_vis.obs["spot_id"] = adata_vis.obs_names
adata_vis.obs["barcode"] = barcodes.reindex(adata_vis.obs_names)["barcode"].values
adata_vis.obs["x"] = barcodes.reindex(adata_vis.obs_names)["x"].values
adata_vis.obs["y"] = barcodes.reindex(adata_vis.obs_names)["y"].values

# 空间坐标
adata_vis.obsm["spatial"] = adata_vis.obs[["x", "y"]].astype(float).values

print("Final spatial AnnData:", adata_vis.shape)
print("\nFirst obs:")
print(adata_vis.obs.head())

print("\nFirst genes:")
print(adata_vis.var_names[:10].tolist())

adata_vis.write(out_path)
print("\nSaved:", out_path)
