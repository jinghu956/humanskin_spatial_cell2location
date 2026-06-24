import scanpy as sc  # 导入 Scanpy，用来读取和处理空间转录组 h5ad 数据。
import pandas as pd  # 导入 pandas，用来读取和保存 csv 格式的 reference signature 和丰度矩阵。
import torch  # 导入 PyTorch，用来检查 GPU/CUDA 并支持 cell2location 模型训练。
from cell2location.models import Cell2location  # 导入 cell2location 的空间定位模型。

torch.set_float32_matmul_precision("high")  # 设置 float32 矩阵乘法精度，提高 GPU 上矩阵计算表现。

sp_path = "cell2location_run/processed/spatial_rna.h5ad"  # 输入：02 脚本生成的空间 RNA AnnData 文件。
signature_path = "cell2location_run/processed/inf_aver_broad_adult.csv"  # 输入：05 脚本生成的细胞类型 reference signature。

model_dir = "cell2location_run/models/cell2location_broad_adult"  # 输出：训练好的 Cell2location 模型保存目录。
out_path = "cell2location_run/results/spatial_cell2location_broad_adult.h5ad"  # 输出：带细胞丰度结果的空间 h5ad 文件。
out_abundance_csv = "cell2location_run/results/cell_abundance_broad_adult.csv"  # 输出：每个 spot 的细胞类型丰度 csv 文件。

adata_vis = sc.read_h5ad(sp_path)  # 读取空间 RNA AnnData。
inf_aver = pd.read_csv(signature_path, index_col=0)  # 读取 reference signature，第一列作为基因名索引。

print("Spatial shape before gene intersection:", adata_vis.shape)  # 打印取共同基因前空间数据维度。
print("Signature shape before gene intersection:", inf_aver.shape)  # 打印取共同基因前 signature 表维度。

# 只保留空间数据和 reference signature 的共同基因。
intersect = adata_vis.var_names.intersection(inf_aver.index)  # 找出空间数据基因和 signature 基因的交集。
adata_vis = adata_vis[:, intersect].copy()  # 空间 AnnData 只保留共同基因，并复制成新对象。
inf_aver = inf_aver.loc[intersect, :].copy()  # signature 表只保留共同基因，并保持和空间数据相同的基因集合。

print("Spatial shape after gene intersection:", adata_vis.shape)  # 打印取共同基因后空间数据维度。
print("Signature shape after gene intersection:", inf_aver.shape)  # 打印取共同基因后 signature 表维度。

print("\nCUDA available:", torch.cuda.is_available())  # 打印 PyTorch 是否检测到可用 CUDA/GPU。
if torch.cuda.is_available():  # 如果当前环境可以使用 GPU。
    print("GPU device count:", torch.cuda.device_count())  # 打印可用 GPU 数量。
    print("Using GPU:", torch.cuda.get_device_name(0))  # 打印第 0 块 GPU 的名称。

Cell2location.setup_anndata(adata=adata_vis)  # 告诉 cell2location 这个 AnnData 是要用于空间模型训练的数据。

# N_cells_per_location 是每个空间点预期包含的细胞数。
# 第一版先用 30，后续可以根据平台和组织厚度调整。
mod = Cell2location(  # 初始化 Cell2location 空间细胞丰度推断模型。
    adata_vis,  # 传入空间 RNA AnnData。
    cell_state_df=inf_aver,  # 传入 reference signature，告诉模型每种细胞类型的基因表达特征。
    N_cells_per_location=30,  # 设置每个空间 spot 预期约包含 30 个细胞。
    detection_alpha=20  # 设置检测效率先验参数，用于建模不同 spot 的测序检测差异。
)

print("\nStart training Cell2location...")  # 打印空间模型训练开始提示。

mod.train(  # 开始训练 Cell2location 模型。
    max_epochs=30000,  # 最多训练 30000 个 epoch，cell2location 通常用较大 epoch 上限并可提前收敛。
    batch_size=None,  # 不分 batch，使用全部 spot 一起训练；这里 spot 数较少可以这样设置。
    train_size=1,  # 使用全部空间 spot 作为训练集，不额外划分验证集。
    accelerator="gpu" if torch.cuda.is_available() else "cpu"  # 有 GPU 就用 GPU，没有 GPU 就用 CPU。
)

print("Training finished.")  # 打印训练完成提示。

print("\nExport posterior...")  # 打印 posterior 导出开始提示。

adata_vis = mod.export_posterior(  # 从训练好的模型导出 posterior 结果并写回空间 AnnData。
    adata_vis,  # 要写入结果的空间 AnnData。
    sample_kwargs={  # posterior 采样相关参数。
        "num_samples": 1000,  # 每个参数抽取 1000 个 posterior samples。
        "batch_size": adata_vis.n_obs,  # 导出时 batch 大小设为空间 spot 总数。
        "accelerator": "gpu" if torch.cuda.is_available() else "cpu"  # posterior 导出同样优先使用 GPU。
    }
)

mod.save(model_dir, overwrite=True)  # 保存训练好的 Cell2location 模型，已有目录则覆盖。
adata_vis.write(out_path)  # 保存包含 posterior 和细胞丰度结果的空间 AnnData。

# 保存每个 spot 的细胞丰度矩阵。
abundance = adata_vis.obsm["q05_cell_abundance_w_sf"]  # 取出 q05 分位数的 spot x cell type 细胞丰度估计矩阵。
abundance.to_csv(out_abundance_csv)  # 把细胞丰度矩阵保存为 csv 文件。

print("\nSaved Cell2location model:", model_dir)  # 打印模型保存目录。
print("Saved result h5ad:", out_path)  # 打印结果 h5ad 保存路径。
print("Saved abundance csv:", out_abundance_csv)  # 打印细胞丰度 csv 保存路径。
print("\nAbundance shape:", abundance.shape)  # 打印丰度矩阵维度：spot 数 x 细胞类型数。
print(abundance.head())  # 打印丰度矩阵前几行，方便快速检查输出格式。
