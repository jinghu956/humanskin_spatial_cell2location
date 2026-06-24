import pandas as pd
import matplotlib.pyplot as plt
import os

infile = "cell2location_run/results/marker_abundance_correlation_summary.csv"
outdir = "cell2location_run/results/figures/correlation_plots"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile)

# 按 relative spearman 排序
df = df.sort_values("spearman_relative_abundance_vs_marker", ascending=True)

plt.figure(figsize=(8, 6))
plt.barh(df["celltype"], df["spearman_relative_abundance_vs_marker"])
plt.xlabel("Spearman correlation\n(relative abundance vs marker score)")
plt.ylabel("Cell type")
plt.title("Marker support for cell2location relative abundance")
plt.tight_layout()

outpath = os.path.join(outdir, "correlation_summary_barplot.png")
plt.savefig(outpath, dpi=300)
plt.close()

print("Saved:", outpath)
print(df[["celltype", "spearman_relative_abundance_vs_marker",
          "pearson_relative_abundance_vs_marker"]])
