#!/bin/bash

# 输入目录
INPUT_DIR="GWAS_AP/05_gwas_mag"
# 输出文件
OUTPUT_FILE="combined_suggestive_snps.txt"

# 定义提示性P值阈值
SUGGESTIVE_P_THRESHOLD=1e-5

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 目录 '$INPUT_DIR' 不存在"
    exit 1
fi

# 获取所有.linear文件
FILES=("$INPUT_DIR"/*.linear)

# 检查是否有文件
if [ ${#FILES[@]} -eq 0 ]; then
    echo "错误: 在 '$INPUT_DIR' 中没有找到.linear文件"
    exit 1
fi

echo "找到 ${#FILES[@]} 个GWAS结果文件"

# 创建输出文件，先写入表头
# 从第一个文件中提取表头，并在开头添加MAG列
echo "提取表头并添加MAG列..."
HEADER=$(head -n 1 "${FILES[0]}")
echo -e "MAG\t$HEADER" > "$OUTPUT_FILE"

# 计数器
COUNTER=0
TOTAL_SNPS=0

# 循环处理每个MAG的GWAS结果
for file in "${FILES[@]}"; do
    # 提取MAG编号
    MAG_NUM=$(basename "$file" | grep -o 'MAG[0-9]\+')
    
    # 检查文件是否为空
    if [ ! -s "$file" ]; then
        echo "跳过空文件: $file"
        continue
    fi
    
    # 获取文件行数（不包括表头）
    FILE_LINES=$(tail -n +2 "$file" | wc -l)
    
    # 提取提示性显著的SNPs (P < 1e-6) 并在每行前添加MAG编号
    # 注意：P值在第15列
    SNPS=$(tail -n +2 "$file" | awk -v p_thresh="$SUGGESTIVE_P_THRESHOLD" -v mag="$MAG_NUM" '
    $15 < p_thresh {
        print mag "\t" $0
    }')
    
    # 计算提取到的SNP数量
    SNP_COUNT=$(echo "$SNPS" | wc -l)
    
    if [ "$SNP_COUNT" -gt 0 ]; then
        echo "从 $MAG_NUM 提取到 $SNP_COUNT 个提示性SNPs (P < $SUGGESTIVE_P_THRESHOLD), 文件共有 $FILE_LINES 个SNPs"
        echo "$SNPS" >> "$OUTPUT_FILE"
        TOTAL_SNPS=$((TOTAL_SNPS + SNP_COUNT))
    else
        echo "$MAG_NUM: 没有找到提示性SNPs (P < $SUGGESTIVE_P_THRESHOLD), 文件共有 $FILE_LINES 个SNPs"
    fi
    
    COUNTER=$((COUNTER + 1))
done

echo ""
echo "========================================"
echo "处理完成!"
echo "已处理 $COUNTER 个MAG文件"
echo "共提取 $TOTAL_SNPS 个提示性显著SNPs (P < $SUGGESTIVE_P_THRESHOLD)"
echo "结果保存在: $OUTPUT_FILE"
echo ""
echo "输出文件信息:"
wc -l "$OUTPUT_FILE"
echo "========================================"

# 可选：如果需要按P值排序
echo "正在按P值排序..."
head -n 1 "$OUTPUT_FILE" > "${OUTPUT_FILE%.txt}_sorted.txt"
tail -n +2 "$OUTPUT_FILE" | sort -k16,16g >> "${OUTPUT_FILE%.txt}_sorted.txt"
echo "已创建排序版本: ${OUTPUT_FILE%.txt}_sorted.txt"
