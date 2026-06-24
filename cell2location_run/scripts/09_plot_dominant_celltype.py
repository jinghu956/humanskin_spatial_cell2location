import pandas as pd
import matplotlib.pyplot as plt
import os

infile = "cell2location_run/results/spot_dominant_celltype.csv"
outdir = "cell2location_run/results/figures"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile, index_col=0)

celltypes = sorted(df["dominant_celltype"].dropna().unique())
color_map = {ct: i for i, ct in enumerate(celltypes)}
colors = df["dominant_celltype"].map(color_map)

plt.figure(figsize=(6, 5))
sc = plt.scatter(
    df["x"],
    df["y"],
    c=colors,
    s=18,
    cmap="tab20"
)
plt.gca().invert_yaxis()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Dominant cell type per spot")

cbar = plt.colorbar(sc, ticks=range(len(celltypes)))
cbar.ax.set_yticklabels(celltypes)

plt.tight_layout()
outpath = os.path.join(outdir, "dominant_celltype_spatial_map.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("Saved:", outpath)

print("\nDominant cell type counts:")
print(df["dominant_celltype"].value_counts())
