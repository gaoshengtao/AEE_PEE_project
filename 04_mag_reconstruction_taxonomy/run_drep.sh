#!/bin/bash
set -euo pipefail

# 输入 MAG 目录
MAG_DIR=/home/gao/data/host_microbe/MAGs/workflow_results/1.MAGs_result/MAGs_assmblyresult/filterMAGs

# 输出目录
OUT_DIR=/home/gao/data/host_microbe/MAGs/drep

# 原始质量信息文件（你已经准备好的）
QUAL_IN=/home/gao/data/host_microbe/MAGs/MAGs_quality.csv

# 生成一个“带后缀”的 genomeInfo 文件，让 dRep 用这个
QUAL_OUT="$OUT_DIR/MAGs_quality_with_suffix.csv"

# 创建输出目录
mkdir -p "$OUT_DIR"

echo "📌 dRep 去冗余开始"
echo "输入 MAG 文件目录: $MAG_DIR"
echo "输出目录: $OUT_DIR"
echo "原始 genomeInfo: $QUAL_IN"
echo "补全后缀的 genomeInfo: $QUAL_OUT"
echo

#############################
# 1. 生成带后缀的 genomeInfo
#############################

if [[ ! -f "$QUAL_IN" ]]; then
    echo "❌ 错误：找不到质量信息文件 $QUAL_IN"
    exit 1
fi

echo "👉 根据实际文件自动补全 genome 列的后缀（.fa/.fna）..."

# 读取原始 CSV，假定为逗号分隔：
# 第一行是表头：genome,completeness,contamination
{
    # 先读表头，原样写出
    IFS= read -r header
    echo "$header"

    # 再逐行处理数据
    while IFS=',' read -r genome completeness contamination rest; do
        # 跳过空行
        [[ -z "$genome" ]] && continue
        # 跳过表头（防止重复）
        if [[ "$genome" == "genome" ]]; then
            continue
        fi

        # 清洗 genome 名：去掉引号、回车、首尾空格
        genome_clean=$(echo "$genome" | tr -d '\r"' | tr -d "'" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')

        # 如果本身已经带后缀且文件存在，就直接用
        if [[ -f "$MAG_DIR/$genome_clean" ]]; then
            genome_file="$genome_clean"
        else
            genome_file=""
            # 尝试补各种可能后缀
            for ext in fa fna fasta; do
                if [[ -f "$MAG_DIR/${genome_clean}.${ext}" ]]; then
                    genome_file="${genome_clean}.${ext}"
                    break
                fi
            done
        fi

        if [[ -z "$genome_file" ]]; then
            echo "⚠ 警告：在 $MAG_DIR 中找不到与 \"$genome_clean\" 对应的 MAG 文件，跳过这一行" >&2
            continue
        fi

        # 把补全后的文件名写入新的 CSV
        if [[ -n "${rest:-}" ]]; then
            echo "$genome_file,$completeness,$contamination,$rest"
        else
            echo "$genome_file,$completeness,$contamination"
        fi

    done
} < "$QUAL_IN" > "$QUAL_OUT"

echo "✅ 补全完成，输出 genomeInfo 文件：$QUAL_OUT"
echo

#####################################
# 2. 收集所有真实存在的 MAG 文件
#####################################

shopt -s nullglob  # 让 *.fa 在没有匹配时展开为空，而不是字面量

genomes=( "$MAG_DIR"/*.fa "$MAG_DIR"/*.fna "$MAG_DIR"/*.fasta )

if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "❌ 错误：在 $MAG_DIR 中没有找到任何 .fa/.fna/.fasta 文件"
    exit 1
fi

echo "检测到的 MAG 文件数量：${#genomes[@]}"
echo "示例前 5 个："
printf '  %s\n' "${genomes[@]:0:5}"
echo
echo "开始运行 dRep..."
echo

#####################################
# 3. 运行 dRep 去冗余（使用质量信息）
#####################################

dRep dereplicate "$OUT_DIR" \
    -g "${genomes[@]}" \
    -p 64 \
    --S_algorithm fastANI \
    --genomeInfo "$QUAL_OUT" \
    --completeness 50 \
    --contamination 10

echo
echo "🎉 dRep 去冗余完成！"
echo "结果已输出到：$OUT_DIR"
echo "代表基因组在：$OUT_DIR/dereplicated_genomes/"

