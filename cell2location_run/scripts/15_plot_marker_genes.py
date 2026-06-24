import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import sparse

adata = sc.read_h5ad("cell2location_run/processed/spatial_rna.h5ad")

outdir = "cell2location_run/results/figures/marker_genes"
os.makedirs(outdir, exist_ok=True)

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

genes_available = set(adata.var_names)

for celltype, genes in marker_genes.items():
    valid_genes = [g for g in genes if g in genes_available]
    print(celltype, "valid markers:", valid_genes)

    for gene in valid_genes:
        idx = adata.var_names.get_loc(gene)
        x = adata.X[:, idx]

        if sparse.issparse(x):
            values = np.asarray(x.toarray()).ravel()
        else:
            values = np.asarray(x).ravel()

        plt.figure(figsize=(5, 5))
        scp = plt.scatter(
            adata.obs["x"],
            adata.obs["y"],
            c=values,
            s=18
        )
        plt.gca().invert_yaxis()
        plt.colorbar(scp, label=f"{gene} count")
        plt.title(f"{gene} ({celltype})")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.tight_layout()

        outpath = os.path.join(outdir, f"{celltype}_{gene}.png")
        plt.savefig(outpath, dpi=300)
        plt.close()

print("Saved marker gene figures to:", outdir)
