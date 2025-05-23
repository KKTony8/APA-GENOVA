### Step 1：转换成loop.bedpe文件
```

#路径设定
INPUT_DIR=/data5/Wuky/CHIFQ/software/YOResult/itx
OUTPUT_DIR=${INPUT_DIR}/loops_bedpe
mkdir -p "$OUTPUT_DIR"

#取前6列，把上游.itx文件之后转换成loop.bedpe文件
for file in "$INPUT_DIR"/*.DA.itx; do
    filename=$(basename "$file")
    sample="${filename%.DA.itx}"
    cut -f1-6 "$file" > "${OUTPUT_DIR}/${sample}.loops.bedpe"
done

```
### Step 2：获取所有loop区域
```

#找出两个年轻和年老bedpe文件的重叠区域，注意顺序
pairToPair -a yHSC-4.loops.bedpe -b oHSC-78.loops.bedpe -type both > yHSC_vs_oHSC_overlap.bedpe

#合并重叠区域，start取小，end取大，保留两个anchor的chr列
awk 'BEGIN{OFS="\t"}{
  chr1 = $1;       # 用第1列作为左端染色体
  chr2 = $4;       # 用第4列作为右端染色体
  start1 = ($2 < $8) ? $2 : $8;
  end1   = ($3 > $9) ? $3 : $9;
  start2 = ($5 < $11) ? $5 : $11;
  end2   = ($6 > $12) ? $6 : $12;
  print chr1, start1, end1, chr2, start2, end2;
}' yHSC_vs_oHSC_overlap.bedpe > merged_overlap_anchors.bedpe
```

```
#取前6列为年轻样本的重叠区域，取后6列为年老样本的重叠区域
cut -f1-6 yHSC_vs_oHSC_overlap.bedpe > overlap_yHSC.bedpe
cut -f7-12 yHSC_vs_oHSC_overlap.bedpe > overlap_oHSC.bedpe

#用两个样本区域分别减去重叠区域得到特有区域
grep -vFf overlap_yHSC.bedpe yHSC-4.loops.bedpe > yHSC-4_specific.bedpe
grep -vFf overlap_oHSC.bedpe oHSC-78.loops.bedpe > oHSC-78_specific.bedpe

#最后cat重叠区域和两个特有区域
cat merged_overlap_anchors.bedpe oHSC-78_specific.bedpe yHSC-4_specific.bedpe > all_loops_merged.bedpe
```

### Step3： 将每个样本的 reads 匹配到 loop 区域上
```
# 设置路径
READS_DIR=/data5/Wuky/CHIFQ/software/YOResult/bam/queryname_sorted/bedpe_files
LOOPS_FILE=/data5/Wuky/CHIFQ/software/YOResult/itx/loops_bedpe/all_loops_merged.bedpe
OUTPUT_DIR=/data5/Wuky/CHIFQ/software/YOResult/read1loop

# 创建输出文件夹
mkdir -p "$OUTPUT_DIR"

# 指定要处理的两个样本
for sample in oHSC-78 yHSC-4; do
    reads_file="${READS_DIR}/${sample}.reads.bedpe"
    output_file="${OUTPUT_DIR}/${sample}.reads_in_loops.txt"

    echo "处理样本：$sample"

    # 使用 bedtools pairtobed 计算每个 read 与 loop 的配对关系
    bedtools pairtobed -a "$reads_file" -b "$LOOPS_FILE" > "$output_file"
done

echo "指定样本的配对处理完成，结果保存在 $OUTPUT_DIR"
```

###Step 4: 构建 count 矩阵
```
# 1. 指定需要合并的两个文件
files=("oHSC-78.reads_in_loops.txt" "yHSC-4.reads_in_loops.txt")

# 2. 设置表头
header="Loop_ID"
samples=""
for file in "${files[@]}"; do
    sample=$(basename "$file" .reads_in_loops.txt)
    header+="\t$sample"
    samples+="$sample "
done
echo -e "$header" > loop_count_matrix.txt

# 3. 统计每个 loop 的 read 数，生成长表
for file in "${files[@]}"; do
    sample=$(basename "$file" .reads_in_loops.txt)
    awk -v SAMPLE="$sample" '
    BEGIN{OFS="\t"}
    {
        loop_id = $11":"$12"-"$13"__"$14":"$15"-"$16
        count[loop_id]++
    }
    END {
        for (l in count) {
            print l, SAMPLE, count[l]
        }
    }' "$file"
done > all_counts.tsv

# 4. 合并为宽表格式矩阵
awk -F'\t' -v SAMPLES="$samples" '
{
    loop=$1; sample=$2; count=$3
    data[loop][sample] = count
    loops[loop]
}
END {
    split(SAMPLES, col_order, " ")
    for (loop in loops) {
        printf "%s", loop
        for (i = 1; i <= length(col_order); i++) {
            s = col_order[i]
            printf "\t%s", (loop in data && s in data[loop]) ? data[loop][s] : 0
        }
        print ""
    }
}' all_counts.tsv >> loop_count_matrix.txt
```
# 5. 删除中间文件
```
rm -f all_counts.tsv

###Step5：均一化使得年轻的coun总和和年老的count的总和一样，计算 log2FC 并筛选差异 loops（|log2FC| > 1）
awk 'BEGIN {
    total1 = 5392163;
    total2 = 4676031;
    pc = 1;  # pseudocount to avoid division by zero
    print "Loop_ID\toHSC-78_CPM\tyHSC-4_CPM\tlog2FC"
}
NR > 1 {
    norm1 = (($2 + pc) / total1) * 1000000;
    norm2 = (($3 + pc) / total2) * 1000000;
    log2fc = log(norm1 / norm2) / log(2);
    if (log2fc > 1 || log2fc < -1)
        printf "%s\t%.4f\t%.4f\t%.4f\n", $1, norm1, norm2, log2fc;
}' loop_count_matrix.txt > loop_count_matrix_log2FC_gt1.txt
```

### Step6：根据 log2FC 分为上调和下调 loops，输出 .bedpe
```
 awk 'NR == 1 { next }
> {
>     split($1, parts, "__");
>     split(parts[1], a, "[:-]");
>     split(parts[2], b, "[:-]");
>     if ($4 >= 1)
>         print a[1], a[2], a[3], b[1], b[2], b[3] >> "old.loop.bedpe";
>     else if ($4 <= -1)
>         print a[1], a[2], a[3], b[1], b[2], b[3] >> "young.loop.bedpe";
> }' loop_count_matrix_log2FC_gt1.txt

```
