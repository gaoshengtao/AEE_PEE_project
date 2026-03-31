#!/usr/bin/env bash
#set -euo pipefail

############################################
# 参数按需调整
############################################

# 已重命名为 .fna 的 MAG 目录（你已有的目录）
MAG_DIR=/home/gao/data/host_microbe/MAGs/gtdbtk_out/genomes_fna

# GTDB-Tk 批处理总输出目录
OUT_DIR=/home/gao/data/host_microbe/MAGs/gtdbtk_out/batch_run

# 每批处理的基因组数量
BATCH_SIZE=500

# 总线程
CPUS=32

# pplacer 用线程（建议 1）
PPLACER_CPUS=1

# scratch 目录
SCRATCH_BASE=/home/gao/data/host_microbe/MAGs/gtdbtk_out/scratch

echo "📌 使用改后缀的 .fna MAGs 进行 GTDB-Tk 分批注释"
echo "MAG_DIR:       $MAG_DIR"
echo "输出目录:      $OUT_DIR"
echo "每批数量:      $BATCH_SIZE"
echo

############################################
# 基本检查
############################################
if ! command -v gtdbtk &>/dev/null; then
    echo "❌ 找不到 gtdbtk，请先激活环境"
    exit 1
fi

if [[ -z "${GTDBTK_DATA_PATH:-}" ]]; then
    echo "❌ GTDBTK_DATA_PATH 未设置"
    exit 1
fi

echo "GTDBTK_DATA_PATH = $GTDBTK_DATA_PATH"
echo

############################################
# 创建输出目录
############################################
mkdir -p "$OUT_DIR"
mkdir -p "$SCRATCH_BASE"

BATCH_LIST_DIR="$OUT_DIR/batch_lists"
ALL_LIST="$BATCH_LIST_DIR/all_genomes.list"

rm -rf "$BATCH_LIST_DIR"
mkdir -p "$BATCH_LIST_DIR"

############################################
# 收集 .fna 文件
############################################
shopt -s nullglob
genomes=( "$MAG_DIR"/*.fna "$MAG_DIR"/*.fa "$MAG_DIR"/*.fasta )

if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "❌ 在 $MAG_DIR 中未找到 fna/fa/fasta 文件"
    exit 1
fi

echo "检测到 MAG 数量：${#genomes[@]}"
echo

############################################
# 生成 genome 名称列表（只保留文件名）
############################################
> "$ALL_LIST"
for f in "${genomes[@]}"; do
    basename "$f" >> "$ALL_LIST"
done

echo "前 5 个 MAG："
head -n 5 "$ALL_LIST"
echo

############################################
# 按行数切分为多个 batch_*
############################################
echo "切分为每批 $BATCH_SIZE 个..."
split -l "$BATCH_SIZE" "$ALL_LIST" "$BATCH_LIST_DIR/batch_"

batch_files=( "$BATCH_LIST_DIR"/batch_* )
echo "共生成批次数：${#batch_files[@]}"
echo

############################################
# 分批运行 GTDB-Tk classify_wf
# 每一批在独立目录创建一个 genomes 子目录，并用软链接指向 MAG_DIR 下的文件
############################################
batch_index=0
for list_file in "${batch_files[@]}"; do
    batch_name=$(basename "$list_file")
    batch_index=$((batch_index + 1))

    batch_out="$OUT_DIR/$batch_name"
    batch_genome_dir="$batch_out/genomes"
    batch_log="$batch_out/run.log"
    batch_scratch="$SCRATCH_BASE/$batch_name"

    mkdir -p "$batch_out"
    mkdir -p "$batch_genome_dir"
    mkdir -p "$batch_scratch"

    # 根据列表创建这一批的软链接
    num=$(wc -l < "$list_file")
    echo "🚀 运行批次 [$batch_index/${#batch_files[@]}] — $batch_name （$num 个基因组）"
    echo "   输出目录:   $batch_out"
    echo "   基因组目录: $batch_genome_dir"
    echo

    # 清理旧链接（防止重跑时残留）
    find "$batch_genome_dir" -type l -delete 2>/dev/null || true

    while read -r fname; do
        src="$MAG_DIR/$fname"
        if [[ ! -f "$src" ]]; then
            echo "⚠️ 警告：找不到文件 $src，跳过"
            continue
        fi
        ln -sf "$src" "$batch_genome_dir/$fname"
        # 如果不想用软链接，可以改成：
        # cp "$src" "$batch_genome_dir/$fname"
    done < "$list_file"

    gtdbtk classify_wf \
        --genome_dir "$batch_genome_dir" \
        --out_dir "$batch_out" \
        --cpus "$CPUS" \
        --pplacer_cpus "$PPLACER_CPUS" \
        --scratch_dir "$batch_scratch" \
        2>&1 | tee "$batch_log"

    EXIT=${PIPESTATUS[0]}
    echo "退出码：$EXIT"
    echo

    if [[ $EXIT -ne 0 ]]; then
        echo "❌ 批次 $batch_name 失败，停止流程"
        exit $EXIT
    fi
done

echo "🎉 所有批次已完成！"
echo

############################################
# 合并 summary 文件
############################################
FINAL_BAC="$OUT_DIR/gtdbtk.bac120.summary.tsv"
FINAL_AR="$OUT_DIR/gtdbtk.ar53.summary.tsv"

echo "📎 合并所有批次 summary..."

# 合并 bac120
first=1
for list_file in "${batch_files[@]}"; do
    batch_name=$(basename "$list_file")
    f="$OUT_DIR/$batch_name/gtdbtk.bac120.summary.tsv"

    if [[ -s "$f" ]]; then
        if (( first )); then
            cat "$f" > "$FINAL_BAC"
            first=0
        else
            tail -n +2 "$f" >> "$FINAL_BAC"
        fi
    fi
done

# 合并 ar53
first=1
for list_file in "${batch_files[@]}"; do
    batch_name=$(basename "$list_file")
    f="$OUT_DIR/$batch_name/gtdbtk.ar53.summary.tsv"

    if [[ -s "$f" ]]; then
        if (( first )); then
            cat "$f" > "$FINAL_AR"
            first=0
        else
            tail -n +2 "$f" >> "$FINAL_AR"
        fi
    fi
done

echo
echo "🎉 GTDB-Tk 批处理与结果合并完成！"
echo "最终结果："
echo "  - $FINAL_BAC"
echo "  - $FINAL_AR（如存在）"

