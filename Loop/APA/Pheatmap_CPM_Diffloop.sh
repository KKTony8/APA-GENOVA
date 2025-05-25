library(pheatmap)

# 1. 读取数据
df <- read.table("loop_count_matrix_log2FC_gt1.txt", header = TRUE)

# 2. 筛选差异 loops
df_filtered <- df[abs(df$log2FC) > 1, ]

# 3. 添加分组标签
df_filtered$Group <- ifelse(df_filtered$log2FC > 1, "old-specific", "young-specific")

# 4. 提取矩阵并以 log2 变换
log2_mat <- log2(as.matrix(df_filtered[, c("yHSC.4_CPM", "oHSC.78_CPM")]))

# 5. 设置行名为 loop ID
rownames(log2_mat) <- df_filtered$Loop_ID

# 6. 构建 annotation 信息
annotation_row <- data.frame(Group = df_filtered$Group)
rownames(annotation_row) <- rownames(log2_mat)

# 7. 设置颜色映射和 breaks
col_fun <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out = 101)  # 根据 log2 值设定范围 [-1, 1]

summary(log2_mat)

library(RColorBrewer)

# Step 1: 对行进行排序
sorted_idx <- order(annotation_row$Group, decreasing = TRUE)  # old 在前
log2_mat_sorted <- log2_mat[sorted_idx, ]
annotation_row_sorted <- annotation_row[sorted_idx, , drop = FALSE]

# Step 2: 设置颜色映射
col_fun <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(0, 8, length.out = 101)  # 根据你的数据分布设置

# Step 3: 设置注释颜色
annotation_colors <- list(
  Group = c(
    "young-specific" = "#4A90E2",  # 蓝色
    "old-specific" = "#D0021B"     # 红色
  )
)

# Step 4: 画热图
pheatmap(
  mat = log2_mat_sorted,
  annotation_row = annotation_row_sorted,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = col_fun,
  breaks = breaks,
  main = "log2(CPM): Young- vs Old-specific Loops (Old on Top)"
)
