###第一步hic转为多分辨率mcool格式文件
#!/bin/bash

# 设置输入输出目录
IN_DIR="/data5/Wuky/CHIFQ/software/YOResult/hic60"
OUT_DIR="${IN_DIR}/mcool"
mkdir -p "$OUT_DIR"

echo "=== Step 1: Converting .hic to .mcool ==="

# 遍历目录下所有 .hic 文件并转换为 .mcool 格式
for infile in "$IN_DIR"/*.hic; do
    basename=$(basename "$infile" .hic)
    outfile="${OUT_DIR}/${basename}.mcool"

    echo "Converting: $basename.hic --> $basename.mcool"
    hic2cool convert "$infile" "$outfile" -r 0
done

####第二步hicrep处理100kb的mcool文件
# 设置工作目录和样本名
DIRIN="/data5/Wuky/CHIFQ/software/YOResult/hic60/mcool"
OUTDIR="/data5/Wuky/CHIFQ/software/YOResult/hic60/scc_output"
mkdir -p "$OUTDIR"  # 如果输出目录不存在，则创建
array=( oHSC-2 oHSC-7 oHSC-8 yHSC-2 yHSC-3 yHSC-4 )
total=$((${#array[@]} - 1))

# 循环样本对组合（两两配对不重复）
for i in $(seq 0 $total); do
    for j in $(seq 0 $total); do
        if [ $j -gt $i ]; then
            FILEIN1="${DIRIN}/${array[$i]}.mcool"
            FILEIN2="${DIRIN}/${array[$j]}.mcool"
            FILEOUT="${OUTDIR}/${array[$i]}_${array[$j]}_scc.txt"

            echo "Running: ${array[$i]} vs ${array[$j]}"

            hicrep "$FILEIN1" "$FILEIN2" "$FILEOUT" \
                --binSize 100000 \
                --h 10 \
                --dBPMax 5000000 \
                --chrNames 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X 
        fi
    done
done
