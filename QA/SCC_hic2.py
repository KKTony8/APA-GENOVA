###可视化SCC矩阵
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# 输入目录
input_dir = "/data5/Wuky/CHIFQ/software/YOResult/hic60/scc_output"

# 获取所有的 scc.txt 文件
files = [f for f in os.listdir(input_dir) if f.endswith("_scc.txt")]

# 提取样本名称（假设样本名称是文件名的前缀）
samples = sorted(set([f.split("_")[0] for f in files]))

# 创建空的 SCC 矩阵
scc_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# 读取每个文件并填充 SCC 矩阵
for f in files:
    path = os.path.join(input_dir, f)
    name = f.replace("_scc.txt", "")
    s1, s2 = name.split("_")

    # 读取该文件中的 SCC 值（假设每个文件中只有一个值）
    with open(path) as fh:
        scores = [float(line.strip()) for line in fh if line.strip() and not line.startswith("#")]
    
    avg_scc = np.mean(scores)

    # 填写矩阵（对称）
    scc_matrix.loc[s1, s2] = avg_scc
    scc_matrix.loc[s2, s1] = avg_scc

# 对角线自己与自己相似性为1
np.fill_diagonal(scc_matrix.values, 1.0)

# 使用 seaborn 绘制热图
plt.figure(figsize=(10, 8))
sns.heatmap(scc_matrix, annot=True, cmap="YlGnBu", linewidths=0.5, cbar_kws={'label': 'SCC'})
plt.title("SCC Matrix Heatmap")
plt.savefig("scc_heatmap.png", dpi=300)
plt.show()
