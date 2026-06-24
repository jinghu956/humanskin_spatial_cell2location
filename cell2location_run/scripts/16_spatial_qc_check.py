import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import sparse

adata = sc.read_h5ad("cell2location_run/processed/spatial_rna.h5ad")

outdir = "cell2location_run/results/figures/spatial_qc"
os.makedirs(outdir, exist_ok=True)

# 识别线粒体基因
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

# 计算 QC 指标
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# 保存 QC 表
qc_cols = [
    "spot_id",
    "barcode",
    "x",
    "y",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt"
]

adata.obs[qc_cols].to_csv("cell2location_run/results/spatial_spot_qc_metrics.csv")

print("Spatial AnnData shape:", adata.shape)

print("\nTotal counts summary:")
print(adata.obs["total_counts"].describe())

print("\nDetected genes per spot summary:")
print(adata.obs["n_genes_by_counts"].describe())

print("\nMitochondrial percentage summary:")
print(adata.obs["pct_counts_mt"].describe())

# 每个基因在多少 spot 表达
X = adata.X
if sparse.issparse(X):
    gene_detected_spots = np.asarray((X > 0).sum(axis=0)).ravel()
else:
    gene_detected_spots = np.asarray((X > 0).sum(axis=0)).ravel()

gene_qc = pd.DataFrame({
    "gene": adata.var_names,
    "detected_spots": gene_detected_spots
}).set_index("gene")

gene_qc.to_csv("cell2location_run/results/spatial_gene_qc_metrics.csv")

print("\nGene detected spots summary:")
print(gene_qc["detected_spots"].describe())

# 画空间 QC 图
for key, title, label in [
    ("total_counts", "Spatial total RNA counts", "Total RNA counts"),
    ("n_genes_by_counts", "Detected genes per spot", "Detected genes"),
    ("pct_counts_mt", "Mitochondrial percentage", "MT percentage")
]:
    plt.figure(figsize=(5, 5))
    scp = plt.scatter(
        adata.obs["x"],
        adata.obs["y"],
        c=adata.obs[key],
        s=18
    )
    plt.gca().invert_yaxis()
    plt.colorbar(scp, label=label)
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{key}_spatial.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

# 画分布直方图
for key, title in [
    ("total_counts", "Total RNA counts per spot"),
    ("n_genes_by_counts", "Detected genes per spot"),
    ("pct_counts_mt", "Mitochondrial percentage per spot")
]:
    plt.figure(figsize=(5, 4))
    plt.hist(adata.obs[key], bins=50)
    plt.xlabel(key)
    plt.ylabel("Number of spots")
    plt.title(title)
    plt.tight_layout()
    outpath = os.path.join(outdir, f"{key}_hist.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

# gene detected spots 分布
plt.figure(figsize=(5, 4))
plt.hist(gene_qc["detected_spots"], bins=50)
plt.xlabel("Detected spots per gene")
plt.ylabel("Number of genes")
plt.title("Gene detection across spots")
plt.tight_layout()
outpath = os.path.join(outdir, "gene_detected_spots_hist.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("\nSaved spot QC table:")
print("cell2location_run/results/spatial_spot_qc_metrics.csv")

print("\nSaved gene QC table:")
print("cell2location_run/results/spatial_gene_qc_metrics.csv")

print("\nSaved QC figures to:")
print(outdir)
