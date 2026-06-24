import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
import os

sp_path = "cell2location_run/processed/spatial_rna.h5ad"
abund_path = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"

out_table = "cell2location_run/results/marker_abundance_correlation_summary.csv"
out_score = "cell2location_run/results/marker_module_scores_with_coords.csv"

adata = sc.read_h5ad(sp_path)
abund = pd.read_csv(abund_path, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in abund.columns if c not in meta_cols]

# log-normalize spatial RNA for marker score
adata_log = adata.copy()
sc.pp.normalize_total(adata_log, target_sum=1e4)
sc.pp.log1p(adata_log)

marker_genes = {
    "Keratinocytes": ["KRT14", "KRT5", "KRT1", "KRT10"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN"],
    "T cell": ["TRAC"],
    "Macrophage": ["C1QA", "C1QB"],
    "DC": ["CD1C", "FCER1A", "CLEC10A"],
    "Vascular endothelium": ["VWF", "KDR"],
    "Lymphatic endothelium": ["PROX1", "PDPN"],
    "Mural cell": ["RGS5", "ACTA2", "MCAM"],
    "Melanocyte": ["DCT"],
    "Mast cell": ["TPSAB1"],
    "Schwann cell": ["S100B", "MPZ", "PLP1", "SOX10"],
    "B cell": ["MS4A1", "CD79A", "CD74"],
}

genes_available = set(adata_log.var_names)

score_df = adata.obs[["spot_id", "barcode", "x", "y"]].copy()
summary_rows = []

for ct, genes in marker_genes.items():
    valid = [g for g in genes if g in genes_available]

    if len(valid) == 0:
        print(ct, "no valid marker genes")
        continue

    gene_idx = [adata_log.var_names.get_loc(g) for g in valid]
    X = adata_log.X[:, gene_idx]

    if sparse.issparse(X):
        values = np.asarray(X.mean(axis=1)).ravel()
    else:
        values = np.asarray(X.mean(axis=1)).ravel()

    score_col = f"{ct}_marker_score"
    score_df[score_col] = values

    if ct in abund.columns:
        a = abund.loc[adata_log.obs_names, ct].values
    else:
        print(ct, "not found in abundance columns")
        continue

    # relative abundance
    total_abund = abund.loc[adata_log.obs_names, celltypes].sum(axis=1).values
    rel_a = a / np.maximum(total_abund, 1e-12)

    spearman_abs = spearmanr(a, values, nan_policy="omit")
    pearson_abs = pearsonr(a, values)

    spearman_rel = spearmanr(rel_a, values, nan_policy="omit")
    pearson_rel = pearsonr(rel_a, values)

    summary_rows.append({
        "celltype": ct,
        "valid_markers": ",".join(valid),
        "n_markers": len(valid),
        "spearman_abs_abundance_vs_marker": spearman_abs.statistic,
        "spearman_abs_pvalue": spearman_abs.pvalue,
        "pearson_abs_abundance_vs_marker": pearson_abs.statistic,
        "pearson_abs_pvalue": pearson_abs.pvalue,
        "spearman_relative_abundance_vs_marker": spearman_rel.statistic,
        "spearman_relative_pvalue": spearman_rel.pvalue,
        "pearson_relative_abundance_vs_marker": pearson_rel.statistic,
        "pearson_relative_pvalue": pearson_rel.pvalue,
    })

summary = pd.DataFrame(summary_rows)
summary = summary.sort_values("spearman_relative_abundance_vs_marker", ascending=False)

summary.to_csv(out_table, index=False)
score_df.to_csv(out_score)

print("Saved:")
print(out_table)
print(out_score)

print("\nCorrelation summary:")
print(summary)
