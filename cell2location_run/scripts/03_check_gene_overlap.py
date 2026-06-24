import scanpy as sc

ref_path = "reference_atlas/adult_skin_atlas.h5ad"
sp_path = "cell2location_run/processed/spatial_rna.h5ad"

adata_ref = sc.read_h5ad(ref_path)
adata_vis = sc.read_h5ad(sp_path)

ref_genes = set(adata_ref.var_names.astype(str))
sp_genes = set(adata_vis.var_names.astype(str))
overlap = sorted(ref_genes & sp_genes)

print("Reference shape:", adata_ref.shape)
print("Spatial shape:", adata_vis.shape)
print("Reference genes:", len(ref_genes))
print("Spatial genes:", len(sp_genes))
print("Overlap genes:", len(overlap))

print("\nFirst 20 overlap genes:")
print(overlap[:20])

if len(overlap) < 3000:
    print("\nWARNING: gene overlap is low. Check gene naming or matrix direction.")
else:
    print("\nGene overlap looks usable.")
