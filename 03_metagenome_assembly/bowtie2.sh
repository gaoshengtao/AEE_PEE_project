#!/bin/bash

# 配置路径
BOWTIE2_PATH="$HOME/apps/bowtie2-2.5.4-linux-x86_64/bowtie2"
BOWTIE2BUILD_PATH="$HOME/apps/bowtie2-2.5.4-linux-x86_64/bowtie2-build"
SCAFFOLDS_DIR="/home/gao/data/host_microbe/assembly/metaspades"
HOSTFREE_DIR="/home/gao/data/host_microbe/hostfree"
UNMAPPING_DIR="/home/gao/data/host_microbe/assembly/unmapping"

# 资源参数
THREADS_PER_SAMPLE=24  # 每个样本使用24线程
TOTAL_THREADS=384      # 服务器总线程数
MAX_PARALLEL=$((TOTAL_THREADS / THREADS_PER_SAMPLE))  # 最大并行数=16

# 获取有效样本列表（仅目录名）
SAMPLE_NAMES=()
for SAMPLE_PATH in "${HOSTFREE_DIR}"/H*; do
  SAMPLE=$(basename "${SAMPLE_PATH}")
  # 检查必要文件是否存在
  if [ -f "${SCAFFOLDS_DIR}/${SAMPLE}/scaffolds.fasta" ] && \
     [ -f "${HOSTFREE_DIR}/${SAMPLE}/${SAMPLE}_paired_1.fastq.gz" ] && \
     [ -f "${HOSTFREE_DIR}/${SAMPLE}/${SAMPLE}_paired_2.fastq.gz" ]; then
    SAMPLE_NAMES+=("${SAMPLE}")
  else
    echo "[WARNING] 缺失文件: ${SAMPLE}" >&2
  fi
done

echo "有效样本数: ${#SAMPLE_NAMES[@]}"

# 并行处理函数
process_sample() {
  local SAMPLE=$1
  local SAMPLE_DIR="${UNMAPPING_DIR}/${SAMPLE}"
  mkdir -p "${SAMPLE_DIR}" || return 1

  # 1. 构建索引
  local INDEX_DIR="${SCAFFOLDS_DIR}/${SAMPLE}/index"
  mkdir -p "${INDEX_DIR}"
  if ! ${BOWTIE2BUILD_PATH} \
    "${SCAFFOLDS_DIR}/${SAMPLE}/scaffolds.fasta" \
    "${INDEX_DIR}/${SAMPLE}_index" > "${SAMPLE_DIR}/index.log" 2>&1; then
    echo "[ERROR] 索引构建失败: ${SAMPLE}" >&2
    return 1
  fi

  # 2. 运行Bowtie2
  local LOG_FILE="${SAMPLE_DIR}/bowtie2.log"
  if ! ${BOWTIE2_PATH} \
    -x "${INDEX_DIR}/${SAMPLE}_index" \
    -1 "${HOSTFREE_DIR}/${SAMPLE}/${SAMPLE}_paired_1.fastq.gz" \
    -2 "${HOSTFREE_DIR}/${SAMPLE}/${SAMPLE}_paired_2.fastq.gz" \
    --threads ${THREADS_PER_SAMPLE} \
    --un-conc-gz "${SAMPLE_DIR}/${SAMPLE}_unmapped_%.fastq.gz" \
    -S /dev/null 2> "${LOG_FILE}"; then
    echo "[ERROR] 比对失败: ${SAMPLE}" >&2
    return 1
  fi

  # 3. 结果统计
  echo "${SAMPLE},$(grep 'reads; of these' ${LOG_FILE} | awk '{print $1}'),$(grep 'overall alignment rate' ${LOG_FILE} | awk '{print $1}')" >> "${UNMAPPING_DIR}/summary.csv"

  # 4. 清理索引
  rm -rf "${INDEX_DIR}"
}

# 导出函数和变量
export -f process_sample
export BOWTIE2_PATH BOWTIE2BUILD_PATH SCAFFOLDS_DIR HOSTFREE_DIR UNMAPPING_DIR THREADS_PER_SAMPLE

# 使用GNU Parallel并行执行
parallel -j ${MAX_PARALLEL} \
         --progress \
         --joblog "${UNMAPPING_DIR}/parallel.log" \
         --retries 2 \
         'process_sample {}' ::: "${SAMPLE_NAMES[@]}"