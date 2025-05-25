library(pheatmap)

# 读取数据
df <- read.table("loop_count_matrix_log2FC_gt1.txt", header = TRUE)

# 筛选差异loops
df_filtered <- df[abs(df$log2FC) > 1, ]
df_filtered$Group <- ifelse(df_filtered$log2FC > 1, "old_loop", "young_loop")

# 分别取两组数据
df_old <- df_filtered[df_filtered$Group == "old_loop", ]
df_young <- df_filtered[df_filtered$Group == "young_loop", ]

# 取CPM并log10转换
mat_old <- log10(as.matrix(df_old[, c("yHSC.4_CPM", "oHSC.78_CPM")]) + 1)
mat_young <- log10(as.matrix(df_young[, c("yHSC.4_CPM", "oHSC.78_CPM")]) + 1)

# 行名设为Loop_ID
rownames(mat_old) <- df_old$Loop_ID
rownames(mat_young) <- df_young$Loop_ID

# 对两个子矩阵分别进行行聚类，获得行顺序
hc_old <- hclust(dist(mat_old))
hc_young <- hclust(dist(mat_young))

order_old <- hc_old$order
order_young <- hc_young$order

# 按聚类顺序重新排序
mat_old_sorted <- mat_old[order_old, ]
mat_young_sorted <- mat_young[order_young, ]

# 合并两个矩阵（行合并）
mat_combined <- rbind(mat_old_sorted, mat_young_sorted)

# 创建合并后的分组注释
group_combined <- c(rep("old_loop", nrow(mat_old_sorted)), rep("young_loop", nrow(mat_young_sorted)))
annotation_row <- data.frame(Group = factor(group_combined, levels = c("old_loop", "young_loop")))
rownames(annotation_row) <- rownames(mat_combined)

# 颜色和断点
col_fun <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(0, 3, length.out = 101)

# 画热图，行不聚类（顺序已是聚类顺序），列不聚类
pheatmap(
  mat = mat_combined,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = annotation_row,
  color = col_fun,
  breaks = breaks,
  main = "log10(CPM + 1) of old and young differential loops"
)
