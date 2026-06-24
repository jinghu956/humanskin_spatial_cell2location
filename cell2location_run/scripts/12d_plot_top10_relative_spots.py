import pandas as pd
import matplotlib.pyplot as plt
import os

infile = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
outdir = "cell2location_run/results/figures/spatial_relative_top10"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in df.columns if c not in meta_cols]

total = df[celltypes].sum(axis=1)
rel = df[celltypes].div(total, axis=0)

top_summary = []

for ct in celltypes:
    threshold = rel[ct].quantile(0.90)
    is_top = rel[ct] >= threshold

    top_summary.append({
        "celltype": ct,
        "top10_threshold": threshold,
        "top_spots": int(is_top.sum())
    })

    plt.figure(figsize=(5, 5))

    plt.scatter(
        df.loc[~is_top, "x"],
        df.loc[~is_top, "y"],
        s=12,
        alpha=0.25
    )

    plt.scatter(
        df.loc[is_top, "x"],
        df.loc[is_top, "y"],
        s=22,
        alpha=0.9
    )

    plt.gca().invert_yaxis()
    plt.title(f"{ct}: top 10% relative abundance spots")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()

    safe_ct = ct.replace(" ", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"top10_relative_{safe_ct}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

top_summary = pd.DataFrame(top_summary)
top_summary.to_csv("cell2location_run/results/top10_relative_spots_summary.csv", index=False)

print("Saved top 10% relative abundance spot figures to:", outdir)
print("\nTop 10% summary:")
print(top_summary)
print("\nSaved table: cell2location_run/results/top10_relative_spots_summary.csv")
