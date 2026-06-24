import pandas as pd
import matplotlib.pyplot as plt
import os

summary_path = "cell2location_run/results/cell_abundance_broad_adult_summary.csv"
outdir = "cell2location_run/results/figures"
os.makedirs(outdir, exist_ok=True)

summary = pd.read_csv(summary_path, index_col=0)
summary = summary.sort_values("sum_abundance", ascending=True)

plt.figure(figsize=(7, 5))
plt.barh(summary.index, summary["sum_abundance"])
plt.xlabel("Total inferred abundance")
plt.ylabel("Cell type")
plt.title("Global inferred cell abundance")
plt.tight_layout()

outpath = os.path.join(outdir, "global_cell_abundance_barplot.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("Saved:", outpath)
print(summary)
