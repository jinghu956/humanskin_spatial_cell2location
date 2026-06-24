import pandas as pd
import matplotlib.pyplot as plt
import os

infile = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
outdir = "cell2location_run/results/figures/spatial_abundance"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in df.columns if c not in meta_cols]

print("Cell types:")
print(celltypes)

for ct in celltypes:
    plt.figure(figsize=(5, 5))
    sc = plt.scatter(
        df["x"],
        df["y"],
        c=df[ct],
        s=18
    )
    plt.gca().invert_yaxis()
    plt.colorbar(sc, label=f"{ct} abundance")
    plt.title(ct)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()

    safe_ct = ct.replace(" ", "_").replace("/", "_")
    outpath = os.path.join(outdir, f"spatial_{safe_ct}.png")
    plt.savefig(outpath, dpi=300)
    plt.close()

print("Saved figures to:", outdir)
