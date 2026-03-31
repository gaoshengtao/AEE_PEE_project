#!/bin/bash

# 配置参数
input_dir="/home/gao/data/hsmilk_cp/rawdata/"
output_base="/home/gao/data/hsmilk_cp/cleandata"
trimmomatic_jar="/home/gao/apps/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_file="/home/gao/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
threads=120

# 创建日志目录
mkdir -p "${output_base}/logs"

# 遍历所有R1文件
find "${input_dir}" -name "*.R1.raw.fastq.gz" | while read r1_file; do
    
    # 构建R2文件路径
    r2_file=$(echo "${r1_file}" | sed 's/\.R1\.raw\.fastq\.gz/.R2.raw.fastq.gz/')
    
    # 从文件名提取样本ID (如TAH01N_1HZ)
    sample_id=$(basename "${r1_file}" | awk -F '.' '{print $1}')
    
    # 创建输出目录
    output_dir="${output_base}/${sample_id}"
    mkdir -p "${output_dir}"

    # 生成输出文件前缀
    base_prefix=$(basename "${r1_file}" .R1.raw.fastq.gz)
    output_prefix="${output_dir}/${base_prefix}"

    # 运行Trimmomatic
    echo "[$(date)] 正在处理样本 ${sample_id}"
    
    java -jar "${trimmomatic_jar}" PE \
        -threads ${threads} \
        -phred33 \
        "${r1_file}" \
        "${r2_file}" \
        "${output_prefix}.R1.paired.fastq.gz" \
        "${output_prefix}.R1.unpaired.fastq.gz" \
        "${output_prefix}.R2.paired.fastq.gz" \
        "${output_prefix}.R2.unpaired.fastq.gz" \
        ILLUMINACLIP:"${adapter_file}:2:30:10:2:True" \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2>&1 | tee "${output_base}/logs/${sample_id}.log"

    echo "[$(date)] 样本 ${sample_id} 处理完成"
done

echo "所有样本处理完成！日志文件保存在：${output_base}/logs"