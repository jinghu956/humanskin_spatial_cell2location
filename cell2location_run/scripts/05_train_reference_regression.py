import scanpy as sc  # 导入 Scanpy，用来读取和处理 h5ad 格式的单细胞参考数据。
import pandas as pd  # 导入 pandas，通常用于表格数据处理；本脚本当前没有直接使用。
import torch  # 导入 PyTorch，用来检测 GPU 并支持 cell2location 的模型训练。
from cell2location.models import RegressionModel  # 导入 cell2location 的参考转录组回归模型。

torch.set_float32_matmul_precision("high")  # 设置 float32 矩阵乘法精度，提高 GPU 上矩阵计算表现。

ref_path = "cell2location_run/processed/ref_broad_adult_downsample.h5ad"  # 输入：04 脚本生成的下采样参考 h5ad。
model_dir = "cell2location_run/models/reference_regression_broad_adult"  # 输出：训练好的 RegressionModel 保存目录。
out_ref_path = "cell2location_run/processed/ref_broad_adult_regression_trained.h5ad"  # 输出：带模型 posterior 结果的参考 h5ad。
out_signature_path = "cell2location_run/processed/inf_aver_broad_adult.csv"  # 输出：每种细胞类型的基因表达 signature 表。

labels_key = "annotation_broad_adult"  # obs 中表示细胞类型标签的列名，模型按这列学习每类细胞的表达特征。

adata_ref = sc.read_h5ad(ref_path)  # 读取下采样后的参考单细胞数据。
adata_ref.var_names_make_unique()  # 确保基因名唯一，避免模型训练或保存时因重复基因名报错。

print("Reference shape:", adata_ref.shape)  # 打印参考数据维度：细胞数 x 基因数。
print("Cell type counts:")  # 打印提示文字，下面输出每类细胞数量。
print(adata_ref.obs[labels_key].value_counts())  # 统计并打印每个细胞类型的细胞数。

print("\nCUDA available:", torch.cuda.is_available())  # 打印当前 PyTorch 是否能使用 CUDA/GPU。
if torch.cuda.is_available():  # 如果检测到可用 GPU。
    print("GPU device count:", torch.cuda.device_count())  # 打印可用 GPU 数量。
    print("Using GPU:", torch.cuda.get_device_name(0))  # 打印第 0 块 GPU 的名称。

RegressionModel.setup_anndata(  # 告诉 cell2location 如何从 AnnData 里读取训练需要的信息。
    adata=adata_ref,  # 传入参考 AnnData 对象。
    labels_key=labels_key  # 指定哪一列 obs 作为细胞类型标签。
)

mod = RegressionModel(adata_ref)  # 用准备好的参考数据初始化 RegressionModel。

print("\nStart training RegressionModel...")  # 打印训练开始提示。

mod.train(  # 开始训练参考转录组回归模型。
    max_epochs=250,  # 最多训练 250 个 epoch。
    batch_size=2500,  # 每个 batch 使用 2500 个细胞。
    train_size=1,  # 使用全部参考细胞作为训练集，不额外划分验证集。
    accelerator="gpu" if torch.cuda.is_available() else "cpu"  # 有 GPU 就用 GPU，没有 GPU 就用 CPU。
)

print("Training finished.")  # 打印训练完成提示。

print("\nExport posterior...")  # 打印 posterior 导出开始提示。

adata_ref = mod.export_posterior(  # 从训练好的模型中导出 posterior 估计结果并写回 AnnData。
    adata_ref,  # 要写入结果的参考 AnnData。
    sample_kwargs={  # posterior 采样相关参数。
        "num_samples": 1000,  # 每个参数抽取 1000 个 posterior samples。
        "batch_size": 2500,  # 导出 posterior 时每个 batch 使用 2500 个细胞。
        "accelerator": "gpu" if torch.cuda.is_available() else "cpu"  # posterior 导出同样优先使用 GPU。
    }
)

mod.save(model_dir, overwrite=True)  # 保存训练好的模型到指定目录，已有目录则覆盖。
adata_ref.write(out_ref_path)  # 保存带 posterior 结果的参考 AnnData。

factor_names = adata_ref.uns["mod"]["factor_names"]  # 从模型结果中取出细胞类型/因子名称。
signature_cols = [f"means_per_cluster_mu_fg_{ct}" for ct in factor_names]  # 拼出每个细胞类型对应的 signature 列名。

inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][signature_cols].copy()  # 从 varm 中提取基因 x 细胞类型的平均表达 signature。
inf_aver.columns = factor_names  # 把 signature 表的列名改成更直观的细胞类型名称。
inf_aver.to_csv(out_signature_path)  # 将细胞类型 signature 保存成 csv，供后续空间定位模型使用。

print("\nSaved reference model:", model_dir)  # 打印模型保存目录。
print("Saved trained reference:", out_ref_path)  # 打印训练后参考 h5ad 的保存路径。
print("Saved cell type signature:", out_signature_path)  # 打印 signature csv 的保存路径。
print("\nSignature shape:", inf_aver.shape)  # 打印 signature 表维度：基因数 x 细胞类型数。
print(inf_aver.head())  # 打印 signature 表前几行，方便快速检查结果格式。
