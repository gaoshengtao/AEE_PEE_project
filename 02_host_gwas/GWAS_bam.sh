#!/bin/bash

# 配置参数
input_dir="/home/gao/backup/host_microbe_metabolism/MJ20240523480-SJ2024112200185/raw_data/"
output_base="/home/gao/data/host_microbe/hostgenome/bwa_out/"
trimmomatic_jar="/home/gao/apps/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_file="/home/gao/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
reference_index="/home/gao/databases/mapping_index/Bos_taurus.ARS-UCD1.3.dna.toplevel.fasta"
threads=128

# 创建日志目录
mkdir -p "${output_base}/logs"
processed_log="${output_base}/logs/processed_samples.log"
touch "${processed_log}"

# 检查BWA索引是否存在
if [ ! -f "${reference_index}.bwt" ]; then
    echo "[$(date)] 为参考基因组创建BWA索引..."
    bwa index "${reference_index}"
fi

# 遍历所有R1文件
find "${input_dir}" -name "*.R1.raw.fastq.gz" | while read r1_file; do
    
    # 构建R2文件路径
    r2_file=$(echo "${r1_file}" | sed 's/\.R1\.raw\.fastq\.gz/.R2.raw.fastq.gz/')
    
    # 从文件名提取实际样本ID (G1202, G1234等)
    sample_id=$(basename "${r1_file}" | awk -F '.' '{print $1}' | awk -F '--' '{print $2}')
    
    # 创建输出目录
    output_dir="${output_base}/${sample_id}"
    mkdir -p "${output_dir}"

    # 生成输出文件前缀
    base_prefix="${sample_id}"
    output_prefix="${output_dir}/${base_prefix}"

    # 检查是否已处理完成
    if [ -f "${output_prefix}.sorted.bam" ] && [ -f "${output_prefix}.sorted.bam.bai" ]; then
        echo "[$(date)] 样本 ${sample_id} 已处理完成（存在 ${output_prefix}.sorted.bam 和 .bai 文件），跳过..."
        echo "${sample_id} $(date) - 已存在结果文件，跳过" >> "${processed_log}"
        continue
    fi

    # 检查Trimmomatic输出文件是否存在
    if [ -f "${output_prefix}.R1.paired.fastq.gz" ] && \
       [ -f "${output_prefix}.R2.paired.fastq.gz" ] && \
       [ -f "${output_prefix}.R1.unpaired.fastq.gz" ] && \
       [ -f "${output_prefix}.R2.unpaired.fastq.gz" ]; then
        echo "[$(date)] 样本 ${sample_id} 已存在Trimmomatic输出文件，跳过质量过滤步骤..."
        echo "${sample_id} $(date) - 已存在Trimmomatic输出文件，跳过质量过滤" >> "${processed_log}"
    else
        # 1. 运行Trimmomatic进行质量过滤
        echo "[$(date)] 正在处理样本 ${sample_id} - 质量过滤..."
        
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
            2>&1 | tee "${output_base}/logs/${sample_id}_trimmomatic.log"

        # 检查Trimmomatic是否成功
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
            echo "[$(date)] 错误：样本 ${sample_id} Trimmomatic 失败！" | tee -a "${output_base}/logs/error.log"
            continue
        fi
    fi

    # 2. 使用BWA进行比对
    echo "[$(date)] 正在处理样本 ${sample_id} - BWA比对..."
    
    bwa mem -t ${threads} \
        "${reference_index}" \
        "${output_prefix}.R1.paired.fastq.gz" \
        "${output_prefix}.R2.paired.fastq.gz" \
        2> "${output_base}/logs/${sample_id}_bwa.log" \
        | samtools view -@ ${threads} -bS - \
        > "${output_prefix}.bam"

    # 检查BWA是否成功
    if [ $? -ne 0 ]; then
        echo "[$(date)] 错误：样本 ${sample_id} BWA比对失败！" | tee -a "${output_base}/logs/error.log"
        continue
    fi

    # 3. 对BAM文件进行排序
    echo "[$(date)] 正在处理样本 ${sample_id} - BAM排序..."
    
    samtools sort -@ ${threads} \
        -o "${output_prefix}.sorted.bam" \
        "${output_prefix}.bam" \
        2> "${output_base}/logs/${sample_id}_sort.log"

    # 检查排序是否成功
    if [ $? -ne 0 ]; then
        echo "[$(date)] 错误：样本 ${sample_id} BAM排序失败！" | tee -a "${output_base}/logs/error.log"
        continue
    fi

    # 删除未排序的BAM文件以节省空间
    rm "${output_prefix}.bam"

    # 4. 为排序后的BAM文件创建索引
    echo "[$(date)] 正在处理样本 ${sample_id} - 创建BAM索引..."
    
    samtools index -@ ${threads} \
        "${output_prefix}.sorted.bam" \
        2> "${output_base}/logs/${sample_id}_index.log"

    # 检查索引是否成功
    if [ $? -ne 0 ]; then
        echo "[$(date)] 错误：样本 ${sample_id} BAM索引创建失败！" | tee -a "${output_base}/logs/error.log"
        continue
    fi

    # 5. 删除中间文件以节省空间（可选）
    echo "[$(date)] 正在处理样本 ${sample_id} - 清理中间文件..."
    rm "${output_prefix}".R*.paired.fastq.gz "${output_prefix}".R*.unpaired.fastq.gz

    # 记录处理完成
    echo "${sample_id} $(date) - 处理完成" >> "${processed_log}"
    echo "[$(date)] 样本 ${sample_id} 处理完成"
done

echo "所有样本处理完成！日志文件保存在：${output_base}/logs"
echo "已处理样本记录在：${processed_log}"