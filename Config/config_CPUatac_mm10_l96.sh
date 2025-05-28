### 10:42 2025/5/28
## /data5/Wuky/CHIFQ/ref/mm10/bwa/mm10.blacklist.bed
## /data5/Wuky/CHIFQ/ref/mm10/bwa/mm10.fa

# !/bin/bash
## Config file for CPU pipeline

## 定义参考基因组的位置（fasta）
export refdir=/data5/Wuky/CHIFQ/ref/

# Define read file suffixes
R1suf="_R1.fastq.gz"
R2suf="_R2.fastq.gz"

# ChIA-PET 参数
extbp=500
selfbp=8000
minMapScore=30

# 基因组版本
genome='mm10'

# Peak calling 参数
organism='mm'  # 'hs' for human, 'mm' for mouse, 'dm' for fly
macsq=0.000001
peakext=150
shiftsize=$((-$peakext / 2))

# ChIA-SIG 参数
PET=3

# Linker 序列
linker=ATAGGCGCTGGTTGCGTTGT

# 资源配置参数
NTHREAD=30
hrs=72
GB=72
seed=12639

# 矩阵分辨率
binsize=100000

# 元信息（非运行必要）
run_type="ChIATAC"
ip_factor="chiatac"
cell_type="mouse"

# SIF 容器存放路径（Juicer 和 CPU 分别）
export addpath=/data5/Wuky/CHIFQ/software

# CPU pipeline 主程序路径
cpuprog="/data5/Wuky/CHIFQ/software/CPU-0.0.1a-r2/source/cpu"

# singularity 路径（某些 sudo 环境中需要手动指定）
export singularityPath="/usr/bin/singularity"
