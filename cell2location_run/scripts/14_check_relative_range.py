import pandas as pd

infile = "cell2location_run/results/cell_abundance_broad_adult_with_coords.csv"
out = "cell2location_run/results/cell_abundance_relative_range_summary.csv"

df = pd.read_csv(infile, index_col=0)

meta_cols = ["spot_id", "barcode", "x", "y"]
celltypes = [c for c in df.columns if c not in meta_cols]

total = df[celltypes].sum(axis=1)
rel = df[celltypes].div(total, axis=0)

summary = rel.describe(
    percentiles=[0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
).T

summary["range"] = summary["max"] - summary["min"]
summary["p95_minus_p05"] = summary["95%"] - summary["5%"]
summary["p99_minus_p01"] = summary["99%"] - summary["1%"]
summary = summary.sort_values("p95_minus_p05", ascending=False)

summary.to_csv(out)

print("Relative abundance range summary:")
print(summary)

print("\nSaved:", out)
