import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

infile = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
outdir = "cell2location_run/results/figures/spatial_relative_abundance_zscore"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in df.columns if c not in meta_cols]

total = df[celltypes].sum(axis=1)
rel = df[celltypes].div(total, axis=0)

z = rel.copy()
for ct in celltypes:
    std = rel[ct].std()
    if std == 0:
        z[ct] = 0
    else:
        z[ct] = (rel[ct] - rel[ct].mean()) / std

z_with_coords = pd.concat([df[meta_cols], z], axis=1)
z_with_coords.to_csv("cell2location_run/results/cell_abundance_broad_adult_relative_zscore_with_coords.csv")

print("Relative abundance z-score summary:")
print(z[celltypes].describe().T)

for ct in celltypes:
    values = z[ct].copy()

    vmin = np.nanpercentile(values, 2)
    vmax = np.nanpercentile(values, 98)

    plt.figure(figsize=(5, 5))
    sc = plt.scatter(
        df["x"],
        df["y"],
        c=values,
        s=18,
        vmin=vmin,
        vmax=vmax
    )
    plt.gca().invert_yaxis()
    plt.colorbar(sc, label=f"{ct} relative z-score")
    plt.title(f"{ct} relative abundance z-score")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()

    safe_ct = ct.replace(" ", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"relative_zscore_{safe_ct}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

print("\nSaved z-score figures to:", outdir)
print("Saved table: cell2location_run/results/cell_abundance_broad_adult_relative_zscore_with_coords.csv")
