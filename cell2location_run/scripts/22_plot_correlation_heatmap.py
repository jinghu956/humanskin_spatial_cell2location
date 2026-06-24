import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

infile = "cell2location_run/results/marker_abundance_correlation_summary.csv"
outdir = "cell2location_run/results/figures/correlation_heatmap"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(infile)

# 选择想画的列
plot_cols = [
    "spearman_abs_abundance_vs_marker",
    "pearson_abs_abundance_vs_marker",
    "spearman_relative_abundance_vs_marker",
    "pearson_relative_abundance_vs_marker"
]

# 更好看的列名
rename_map = {
    "spearman_abs_abundance_vs_marker": "Spearman\nAbs",
    "pearson_abs_abundance_vs_marker": "Pearson\nAbs",
    "spearman_relative_abundance_vs_marker": "Spearman\nRelative",
    "pearson_relative_abundance_vs_marker": "Pearson\nRelative"
}

# 按 relative spearman 从高到低排序
df = df.sort_values("spearman_relative_abundance_vs_marker", ascending=False).reset_index(drop=True)

heat = df.set_index("celltype")[plot_cols].copy()
heat = heat.rename(columns=rename_map)

# 保存排序后的表，便于后续查看
heat.to_csv("cell2location_run/results/marker_abundance_correlation_heatmap_values.csv")

data = heat.values
row_labels = heat.index.tolist()
col_labels = heat.columns.tolist()

fig_w = 7
fig_h = max(4, 0.5 * len(row_labels) + 1.5)

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(data, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)

# 坐标轴标签
ax.set_xticks(np.arange(len(col_labels)))
ax.set_yticks(np.arange(len(row_labels)))
ax.set_xticklabels(col_labels, rotation=45, ha="right")
ax.set_yticklabels(row_labels)

# 在每个格子里写数值
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        val = data[i, j]
        text_color = "white" if abs(val) > 0.45 else "black"
        ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color, fontsize=9)

# 标题
ax.set_title("Correlation between cell2location abundance and marker scores")

# colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Correlation coefficient")

plt.tight_layout()

outpath_png = os.path.join(outdir, "correlation_heatmap.png")
outpath_pdf = os.path.join(outdir, "correlation_heatmap.pdf")

plt.savefig(outpath_png, dpi=300)
plt.savefig(outpath_pdf)
plt.close()

print("Saved:")
print(outpath_png)
print(outpath_pdf)
print("cell2location_run/results/marker_abundance_correlation_heatmap_values.csv")

print("\nHeatmap values:")
print(heat)
