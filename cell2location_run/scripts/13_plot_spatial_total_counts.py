import scanpy as sc
import matplotlib.pyplot as plt
import os

adata = sc.read_h5ad("cell2location_run/processed/spatial_rna.h5ad")

outdir = "cell2location_run/results/figures/qc"
os.makedirs(outdir, exist_ok=True)

adata.obs["total_counts"] = adata.X.sum(axis=1)
if not isinstance(adata.obs["total_counts"].iloc[0], float):
    adata.obs["total_counts"] = adata.obs["total_counts"].astype(float)

print("Spatial total counts summary:")
print(adata.obs["total_counts"].describe())

plt.figure(figsize=(5, 5))
sc = plt.scatter(
    adata.obs["x"],
    adata.obs["y"],
    c=adata.obs["total_counts"],
    s=18
)
plt.gca().invert_yaxis()
plt.colorbar(sc, label="Total RNA counts")
plt.title("Spatial total RNA counts per spot")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()

outpath = os.path.join(outdir, "spatial_total_counts.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("Saved:", outpath)
