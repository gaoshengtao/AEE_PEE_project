#!/bin/bash
set -euo pipefail  # 严格错误检查

# --------------------------
# 参数配置区（根据需求修改）
# --------------------------
INPUT_BASE="/home/gao/data/host_microbe/assembly/metaspades"  # 输入根目录
CLUSTERED_DIR="${INPUT_BASE}/clustered_contigs"              # 结果输出根目录（按样本分文件夹）
SIMILARITY=0.90                                              # 去冗余相似性阈值（对应输出后缀）
THREADS=64                                                 # 线程数
MIN_LENGTH=500                                               # 过滤的最小contig长度（bp）

# 定义完成标志文件的后缀（修正：补充下划线，与cd-hit-est生成的.clstr文件一致）
COMPLETION_FLAG_SUFFIX="_group_${SIMILARITY}.fasta.clustered.clstr"

# --------------------------
# 初始化检查与准备
# --------------------------
# 检查输入文件是否存在
contig_files=("${INPUT_BASE}"/*/contigs.fasta)
if [[ ${#contig_files[@]} -eq 0 ]]; then
    echo "错误：未找到任何contigs.fasta文件！"
    exit 1
fi
echo "检测到 ${#contig_files[@]} 个样本的contigs.fasta文件，开始处理..."

# 创建结果根目录（若不存在）
mkdir -p "${CLUSTERED_DIR}"

# 检查seqtk是否安装（关键依赖）
if ! command -v seqtk &> /dev/null; then
    echo "错误：未找到seqtk工具！请先安装seqtk（https://github.com/lh3/seqtk）"
    exit 1
fi

# 检查cd-hit-est是否安装（关键依赖）
if ! command -v cd-hit-est &> /dev/null; then
    echo "错误：未找到cd-hit-est工具！请先安装CD-HIT（http://weizhongli-lab.org/cd-hit/）"
    exit 1
fi

# --------------------------
# 单样本处理核心逻辑（逐个样本过滤+去冗余）
# --------------------------
for contig_file in "${contig_files[@]}"; do
    # ----------------------
    # 步骤1：提取样本信息
    # ----------------------
    sample_dir=$(dirname "${contig_file}")       # 样本所在目录（如 /path/to/H0619）
    sample_name=$(basename "${sample_dir}")      # 样本名（如 H0619）
    echo "===== 处理样本：${sample_name} ====="

    # ----------------------
    # 步骤2：创建样本专属输出目录
    # ----------------------
    sample_output_dir="${CLUSTERED_DIR}/${sample_name}"
    mkdir -p "${sample_output_dir}"
    echo "样本${sample_name}的输出目录：${sample_output_dir}"

    # ----------------------
    # 步骤3：定义完成标志文件路径（修正后路径）
    # ----------------------
    completion_flag="${sample_output_dir}/${sample_name}${COMPLETION_FLAG_SUFFIX}"
    
    # 检查是否已完成分析（标志文件存在则跳过）
    if [[ -f "${completion_flag}" ]]; then
        echo "样本${sample_name}已完成分析（标志文件${completion_flag}存在），跳过本次处理。"
        continue
    fi

    # ----------------------
    # 步骤4：用seqtk过滤长度≥500bp的contigs
    # ----------------------
    filtered_file="${sample_output_dir}/filtered_${MIN_LENGTH}.fasta"
    echo "过滤长度≥${MIN_LENGTH}bp的contigs到：${filtered_file}..."
    seqtk seq -L "${MIN_LENGTH}" "${contig_file}" > "${filtered_file}"

    # 验证过滤后文件是否非空（避免空文件导致cd-hit-est报错）
    if [[ ! -s "${filtered_file}" ]]; then
        echo "警告：样本${sample_name}过滤后文件为空！跳过去冗余步骤。"
        rm -f "${filtered_file}"
        continue
    fi

    # ----------------------
    # 步骤5：用cd-hit-est去冗余（保留90%相似性）
    # ----------------------
    output_file="${sample_output_dir}/${sample_name}_group_${SIMILARITY}.fasta"
    echo "运行cd-hit-est去冗余（相似性${SIMILARITY}）..."
    cd-hit-est \
        -i "${filtered_file}" \
        -o "${output_file}" \
        -c "${SIMILARITY}" \
        -M 0 \
        -T "${THREADS}"

    # ----------------------
    # 步骤6：验证去冗余结果是否生成（确保流程完整性）
    # ----------------------
    if [[ ! -f "${output_file}" ]]; then
        echo "错误：样本${sample_name}去冗余失败，未生成结果文件！"
        rm -f "${filtered_file}"  # 清理临时文件
        continue
    fi

    # ----------------------
    # 步骤7：清理临时过滤文件
    # ----------------------
    rm -f "${filtered_file}"
    echo "已清理临时文件：${filtered_file}"

    # ----------------------
    # 步骤8：创建完成标志文件（标记本次处理成功）
    # ----------------------
    touch "${completion_flag}"
    echo "样本${sample_name}处理完成，已创建完成标志文件：${completion_flag}"
done

echo "===== 所有样本处理完成 ====="
echo "最终结果保存在：${CLUSTERED_DIR}（按样本分文件夹）"