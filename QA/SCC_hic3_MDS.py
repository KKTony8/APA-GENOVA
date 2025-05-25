###根据前面创建的SCC矩阵，绘制MDS降维图
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS

# 输入目录
input_dir = "/data5/Wuky/CHIFQ/software/YOResult/hic60/scc_output"

# 获取所有的 scc.txt 文件
files = [f for f in os.listdir(input_dir) if f.endswith("_scc.txt")]

# 提取样本名称（假设样本名称是文件名前缀）
samples = sorted(set([f.split("_")[0] for f in files] + [f.split("_")[1] for f in files]))

# 创建空的 SCC 矩阵
scc_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# 读取每个文件并填充 SCC 矩阵
for f in files:
    path = os.path.join(input_dir, f)
    name = f.replace("_scc.txt", "")
    s1, s2 = name.split("_")

    with open(path) as fh:
        scores = [float(line.strip()) for line in fh if line.strip() and not line.startswith("#")]
    avg_scc = np.mean(scores)

    scc_matrix.loc[s1, s2] = avg_scc
    scc_matrix.loc[s2, s1] = avg_scc

# 自身相似性设为 1
np.fill_diagonal(scc_matrix.values, 1.0)

# 将 SCC 相似性转换为距离（1 - SCC）
distance_matrix = 1 - scc_matrix

# MDS 分析
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
coords = mds.fit_transform(distance_matrix)

# 画 MDS 散点图
plt.figure(figsize=(8, 6))
for i, sample in enumerate(scc_matrix.index):
    plt.scatter(coords[i, 0], coords[i, 1], s=100)
    plt.text(coords[i, 0] + 0.005, coords[i, 1] + 0.005, sample, fontsize=10)

plt.title("MDS of Samples Based on SCC")
plt.xlabel("MDS1")
plt.ylabel("MDS2")
plt.tight_layout()
plt.savefig("scc_mds.png", dpi=300)
plt.show()
