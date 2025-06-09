####22:55 2025/6/8 wky

###第一步： 转换hic文件：对hic2cool转换.hic文件为.mcool文件
###第二步：查找TAD，运用hicexplorer软件中hicFindTADs

#!/bin/bash

# 创建输出目录（如果不存在）
mkdir -p /data5/Wuky/CHIFQ/software/YOResult/hic60/TADnew/

# 年轻样本 TAD 识别
hicFindTADs \
  --matrix /data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/yHSC-A.mcool::resolutions/50000 \
  --outPrefix /data5/Wuky/CHIFQ/software/YOResult/hic60/TADnew/yHSC_TADs \
  --correctForMultipleTesting fdr \
  --minDepth 300000 \
  --maxDepth 600000

# 年老样本 TAD 识别
hicFindTADs \
  --matrix /data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/oHSC-8.mcool::resolutions/50000 \
  --outPrefix /data5/Wuky/CHIFQ/software/YOResult/hic60/TADnew/oHSC_TADs \
  --correctForMultipleTesting fdr \
  --minDepth 300000 \
  --maxDepth 600000

###第三步： 绘制TAD长度与数量的关系图

library(tidyverse)

# 定义函数读取并处理 TAD 文件
read_tad_file <- function(file_path, sample_name) {
  read_delim(file_path,
             delim = "\t", col_names = FALSE, trim_ws = TRUE) %>%
    select(chr = X1, start = X2, end = X3, id = X4) %>%
    mutate(length = end - start,
           sample = sample_name)
}

# 读取年轻和年老样本的 TAD 数据
young_tads <- read_tad_file("yHSC_TADs_domains.bed", "young")
old_tads   <- read_tad_file("oHSC_TADs_domains.bed", "old")

# 合并两个样本的数据
all_tads <- bind_rows(young_tads, old_tads)

# 计算统计指标
tads_summary <- all_tads %>%
  group_by(sample) %>%
  summarise(
    count = n(),
    max_length = max(length),
    min_length = min(length),
    mean_length = mean(length),
    median_length = median(length),
    total_length = sum(length),
    .groups = "drop"
  )

# 输出结果
print(tads_summary)

####第四步： 绘制TAD长度分布密度图
#密度图
library(tidyverse)

# 转换长度为 Mb
all_tads <- all_tads %>%
  mutate(length_Mb = length / 1e6)

# 密度线图（去掉填充）
ggplot(all_tads, aes(x = length_Mb, color = sample)) +
  geom_density(size = 1.2) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(title = "TAD Length Distribution",
       x = "TAD Length (Mb)",
       y = "Density",
       color = "Sample") +
  theme_minimal()

###*直方图
ggplot(all_tads, aes(x = length / 1e6, fill = sample)) +
  geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.2) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(title = "TAD Length Histogram",
       x = "TAD Length (Mb)",
       y = "Count",
       fill = "Sample") +
  theme_minimal()
