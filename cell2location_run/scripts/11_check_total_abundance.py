import pandas as pd
import matplotlib.pyplot as plt
import os

infile = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
outdir = "cell2location_run/results/figures/qc"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in df.columns if c not in meta_cols]

df["total_abundance"] = df[celltypes].sum(axis=1)

print("Total abundance summary:")
print(df["total_abundance"].describe())

plt.figure(figsize=(5, 5))
sc = plt.scatter(df["x"], df["y"], c=df["total_abundance"], s=18)
plt.gca().invert_yaxis()
plt.colorbar(sc, label="Total inferred abundance")
plt.title("Total inferred cell abundance per spot")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
outpath = os.path.join(outdir, "total_abundance_spatial.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("Saved:", outpath)
