import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import sparse

adata = sc.read_h5ad("cell2location_run/processed/spatial_rna.h5ad")

base_outdir = "cell2location_run/results/figures/marker_genes_relative"
os.makedirs(base_outdir, exist_ok=True)
os.makedirs(os.path.join(base_outdir, "fraction"), exist_ok=True)
os.makedirs(os.path.join(base_outdir, "lognorm"), exist_ok=True)
os.makedirs(os.path.join(base_outdir, "zscore"), exist_ok=True)

marker_genes = {
    "Keratinocytes": ["KRT14", "KRT5", "KRT1", "KRT10"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "T_cell": ["CD3D", "CD3E", "TRAC"],
    "Macrophage": ["LYZ", "C1QA", "C1QB"],
    "DC": ["CD1C", "FCER1A", "CLEC10A"],
    "Vascular_endothelium": ["PECAM1", "VWF", "KDR"],
    "Lymphatic_endothelium": ["PROX1", "LYVE1", "PDPN"],
    "Mural_cell": ["RGS5", "ACTA2", "MCAM"],
    "Melanocyte": ["PMEL", "MLANA", "TYR", "DCT"],
    "Mast_cell": ["TPSAB1", "TPSB2", "CPA3"],
    "Schwann_cell": ["S100B", "MPZ", "PLP1", "SOX10"],
    "B_cell": ["MS4A1", "CD79A", "CD74"]
}

# 原始 total counts
if sparse.issparse(adata.X):
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
else:
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

# 做一个 copy 用于 log-normalization
adata_log = adata.copy()
sc.pp.normalize_total(adata_log, target_sum=1e4)
sc.pp.log1p(adata_log)

genes_available = set(adata.var_names)

for celltype, genes in marker_genes.items():
    valid_genes = [g for g in genes if g in genes_available]
    print(celltype, "valid markers:", valid_genes)

    for gene in valid_genes:
        idx = adata.var_names.get_loc(gene)

        # raw counts
        x_raw = adata.X[:, idx]
        if sparse.issparse(x_raw):
            raw_values = np.asarray(x_raw.toarray()).ravel()
        else:
            raw_values = np.asarray(x_raw).ravel()

        # 1) fraction of total counts
        frac_values = raw_values / np.maximum(total_counts, 1)

        # 2) log-normalized expression
        x_log = adata_log.X[:, idx]
        if sparse.issparse(x_log):
            log_values = np.asarray(x_log.toarray()).ravel()
        else:
            log_values = np.asarray(x_log).ravel()

        # 3) z-score across spots (based on lognorm)
        std = np.std(log_values)
        if std == 0:
            z_values = np.zeros_like(log_values)
        else:
            z_values = (log_values - np.mean(log_values)) / std

        # 画图函数
        def save_plot(values, outdir, suffix, cbar_label, title):
            vmin = np.nanpercentile(values, 2)
            vmax = np.nanpercentile(values, 98)
            if vmin == vmax:
                vmin = np.min(values)
                vmax = np.max(values)

            plt.figure(figsize=(5, 5))
            scp = plt.scatter(
                adata.obs["x"],
                adata.obs["y"],
                c=values,
                s=18,
                vmin=vmin,
                vmax=vmax
            )
            plt.gca().invert_yaxis()
            plt.colorbar(scp, label=cbar_label)
            plt.title(title)
            plt.xlabel("x")
            plt.ylabel("y")
            plt.tight_layout()

            outpath = os.path.join(outdir, f"{celltype}_{gene}_{suffix}.png")
            plt.savefig(outpath, dpi=300)
            plt.close()

        save_plot(
            frac_values,
            os.path.join(base_outdir, "fraction"),
            "fraction",
            f"{gene} / total_counts",
            f"{gene} ({celltype}) fraction"
        )

        save_plot(
            log_values,
            os.path.join(base_outdir, "lognorm"),
            "lognorm",
            f"{gene} log-normalized",
            f"{gene} ({celltype}) log-normalized"
        )

        save_plot(
            z_values,
            os.path.join(base_outdir, "zscore"),
            "zscore",
            f"{gene} z-score",
            f"{gene} ({celltype}) z-score"
        )

print("Saved marker gene relative figures to:", base_outdir)
