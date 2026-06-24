import scanpy as sc  # 导入 Scanpy，用来读取和处理 h5ad 单细胞数据。
import anndata as ad  # 导入 AnnData 类，用来重新构建参考单细胞对象。
import numpy as np  # 导入 NumPy，用来做随机抽样和数值检查。
from scipy import sparse  # 导入稀疏矩阵工具，用来判断表达矩阵是不是 sparse matrix。

ref_path = "reference_atlas/adult_skin_atlas.h5ad"  # 输入参考图谱 h5ad 文件路径。
out_path = "cell2location_run/processed/ref_broad_adult_downsample.h5ad"  # 输出下采样后的参考 h5ad 文件路径。

labels_key = "annotation_broad_adult"  # obs 里保存粗粒度细胞类型注释的列名。
max_cells_per_type = 5000  # 每个细胞类型最多保留的细胞数，超过这个数就随机下采样。
random_state = 0  # 随机种子，保证每次下采样结果可重复。

adata0 = sc.read_h5ad(ref_path)  # 读取原始成人皮肤参考图谱。

print("Original adata:", adata0.shape)  # 打印原始 AnnData 的维度：细胞数 x 基因数。
print("raw exists:", adata0.raw is not None)  # 检查这个对象里是否有 raw 层。

if adata0.raw is None:  # 如果没有 raw 层，说明无法从 raw.X 取原始 counts。
    raise ValueError("adata.raw is None. Cannot find raw counts.")  # 主动报错并停止，避免后面用错表达矩阵。

# 用 raw.X 重建 AnnData，确保 X 是原始 count。
adata_ref = ad.AnnData(  # 新建一个 AnnData，专门作为 cell2location 的参考数据。
    X=adata0.raw.X.copy(),  # 把 raw.X 复制到新的 X，作为原始 count 表达矩阵。
    obs=adata0.obs.copy(),  # 复制细胞元数据，保留每个细胞的注释信息。
    var=adata0.raw.var.copy()  # 复制 raw 里的基因元数据，使其和 raw.X 的基因维度对应。
)

adata_ref.var_names = adata0.raw.var_names.astype(str)  # 把 raw 里的基因名转成字符串并赋给新对象。
adata_ref.var_names_make_unique()  # 如果有重复基因名，自动加后缀使基因名唯一。

# 把 count 转成整数。
if sparse.issparse(adata_ref.X):  # 判断表达矩阵是否是稀疏矩阵。
    adata_ref.X.data = np.rint(adata_ref.X.data).astype("int32")  # 稀疏矩阵只处理非零值：四舍五入后转成 int32。
else:  # 如果表达矩阵是普通 dense 数组。
    adata_ref.X = np.rint(adata_ref.X).astype("int32")  # 对整个矩阵四舍五入并转成 int32。

# 保留有细胞类型注释的细胞。
adata_ref = adata_ref[adata_ref.obs[labels_key].notna()].copy()  # 过滤掉 annotation_broad_adult 缺失的细胞。
adata_ref.obs[labels_key] = adata_ref.obs[labels_key].astype("category")  # 把细胞类型列转为 category，方便后面按类别遍历。

print("Reference using raw counts:", adata_ref.shape)  # 打印过滤后的参考数据维度。
print("\nOriginal cell type counts:")  # 打印提示文字，下面是原始各细胞类型数量。
print(adata_ref.obs[labels_key].value_counts())  # 统计并打印每个细胞类型的细胞数。

# 按细胞类型抽样。
rng = np.random.default_rng(random_state)  # 创建 NumPy 随机数生成器，并固定随机种子。
selected = []  # 用来保存最终选中的细胞下标。

for ct in adata_ref.obs[labels_key].cat.categories:  # 逐个遍历所有细胞类型类别。
    idx = np.where(adata_ref.obs[labels_key].values == ct)[0]  # 找到当前细胞类型对应的所有细胞下标。
    if len(idx) > max_cells_per_type:  # 如果该细胞类型细胞数超过设定上限。
        idx = rng.choice(idx, max_cells_per_type, replace=False)  # 不放回随机抽取最多 5000 个细胞。
    selected.extend(idx)  # 把当前细胞类型保留的下标加入总列表。

adata_ref_ds = adata_ref[selected].copy()  # 根据选中的下标生成下采样后的参考 AnnData。

print("\nDownsampled reference:", adata_ref_ds.shape)  # 打印下采样后的数据维度。
print("\nDownsampled cell type counts:")  # 打印提示文字，下面是下采样后各细胞类型数量。
print(adata_ref_ds.obs[labels_key].value_counts())  # 统计并打印下采样后每个细胞类型的细胞数。

# 再检查一次是否为整数 count。
if sparse.issparse(adata_ref_ds.X):  # 如果下采样后的表达矩阵是稀疏矩阵。
    vals = adata_ref_ds.X.data[:1000000]  # 取最多前 100 万个非零表达值用于检查。
else:  # 如果下采样后的表达矩阵是 dense 数组。
    vals = np.asarray(adata_ref_ds.X).ravel()[:1000000]  # 展平成一维数组后取最多前 100 万个值用于检查。

print("\nCount check:")  # 打印 count 检查的标题。
print("min:", vals.min())  # 打印抽样检查值中的最小值。
print("max:", vals.max())  # 打印抽样检查值中的最大值。
print("first nonzero:", vals[vals != 0][:20])  # 打印前 20 个非零表达值，方便肉眼确认。
print("all sampled values are integers:", np.all(np.equal(vals, np.floor(vals))))  # 检查抽样到的值是否全部为整数。

adata_ref_ds.write(out_path)  # 把下采样后的参考数据写入 h5ad 文件。
print("\nSaved:", out_path)  # 打印保存路径，提示脚本完成。
