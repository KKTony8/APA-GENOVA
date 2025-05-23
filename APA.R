library(GENOVA)

# 路径设定
shared_loop_bedpe <- "/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/old.loop.bedpe"
young_mcool_path  <- "/data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/yHSC-4.mcool"
old_mcool_path    <- "/data5/Wuky/CHIFQ/software/YOResult/hic60/mcool/oHSC-78.mcool"

# 读取 mcool 文件
yHSC_4_10kb <- load_contacts(signal_path = young_mcool_path,
                             sample_name = "yHSC_4",
                             resolution = 10000,
                             balancing = TRUE,
                             colour = "black")

oHSC_78_10kb <- load_contacts(signal_path = old_mcool_path,
                             sample_name = "oHSC_78",
                             resolution = 10000,
                             balancing = TRUE,
                             colour = "red")

# 正确读取并格式化共享 loop 文件（自动识别空格/tab 分隔）
shared_loops <- read.table(shared_loop_bedpe, header=FALSE, sep="", stringsAsFactors=FALSE)
colnames(shared_loops) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

# 去除"chr"前缀（GENOVA使用不带chr的染色体命名）
shared_loops$chr1 <- gsub("^chr", "", shared_loops$chr1)
shared_loops$chr2 <- gsub("^chr", "", shared_loops$chr2)
# 分别对两个样本做 APA 分析
apa_young <- APA(list(yHSC_4 = yHSC_4_10kb),
                 bedpe = shared_loops,
                 dist_thres = c(200e3, Inf))

apa_old <- APA(list(oHSC_78 = oHSC_78_10kb),
               bedpe = shared_loops,
               dist_thres = c(200e3, Inf))
library(ggplot2)
library(dplyr)

# 计算 log2 value
p_old <- visualise(apa_old, title = "APA Heatmap - Old", raw = TRUE)
p_young <- visualise(apa_young, title = "APA Heatmap - Young", raw = TRUE)

p_old$data <- p_old$data %>%
  mutate(value = log2(value + 1))
p_young$data <- p_young$data %>%
  mutate(value = log2(value + 1))

# 获取两个图中 log2(value+1) 的最大值
max_val <- max(c(p_old$data$value, p_young$data$value), na.rm = TRUE)

# 画 Old 图，统一颜色映射范围
p_old + aes(fill = value) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, max_val)) +
  ggtitle("APA Heatmap - Old-Oldloop (log2 scaled)") +
  theme_minimal()

# 画 Young 图，统一颜色映射范围
p_young + aes(fill = value) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, max_val)) +
  ggtitle("APA Heatmap - Young-Oldloop (log2 scaled)") +
  theme_minimal()

# 做差异分析
apa_diff <- APA(list("young"=yHSC_4_10kb, "old"=oHSC_78_10kb),
                bedpe=shared_loops,
                dist_thres=c(200e3, Inf))

# 保存差异的heatmap
visualise(apa_diff,
          title="Young vs Old HSC Loops APA",
          colour_lim=c(0, 40),
          colour_lim_contrast=c(-5, 5),
          metric="diff",
          contrast=1)

#计算APA score
quant <- quantify(apa_young)
apa_score_young <- quant$per_sample$foldchange
cat("APA score (young):", apa_score_young, "\n")

quant_old <- quantify(apa_old)
apa_score_old <- quant_old$per_sample$foldchange
cat("APA score (old):", apa_score_old, "\n")
