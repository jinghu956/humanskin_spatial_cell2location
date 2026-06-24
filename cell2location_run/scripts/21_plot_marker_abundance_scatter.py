import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import spearmanr, pearsonr

abund_file = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
score_file = "cell2location_run/results/marker_module_scores_with_coords.csv"

outdir = "cell2location_run/results/figures/correlation_scatter"
os.makedirs(outdir, exist_ok=True)

abund = pd.read_csv(abund_file, index_col=0)
score = pd.read_csv(score_file, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in abund.columns if c not in meta_cols]

# 相对 abundance
total_abund = abund[celltypes].sum(axis=1)
rel_abund = abund[celltypes].div(total_abund, axis=0)

# 找每个 celltype 对应 marker score 列
for ct in celltypes:
    score_col = f"{ct}_marker_score"
    if score_col not in score.columns:
        print(f"Skip {ct}: no marker score column")
        continue

    x = rel_abund[ct]
    y = score[score_col]

    sp = spearmanr(x, y, nan_policy="omit")
    pr = pearsonr(x, y)

    plt.figure(figsize=(5, 5))
    plt.scatter(x, y, s=10, alpha=0.5)
    plt.xlabel(f"{ct} relative abundance")
    plt.ylabel(f"{ct} marker score")
    plt.title(
        f"{ct}\n"
        f"Spearman={sp.statistic:.3f}, p={sp.pvalue:.2e}\n"
        f"Pearson={pr.statistic:.3f}, p={pr.pvalue:.2e}"
    )
    plt.tight_layout()

    safe_ct = ct.replace(" ", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"{safe_ct}_relative_abundance_vs_marker_score.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

    print("Saved:", outpath)
