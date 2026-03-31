#!/bin/bash

# 修复版 (已验证)
input_base="/home/gao/data/host_microbe/hostfree/"
output_dir="/home/gao/data/host_microbe/assembly/megahit"
threads=360
batch_size=50
max_parallel=2

# 关键配置1：指定真实路径（根据实际路径修改）
MEGAHIT_PATH="/home/gao/miniconda3/envs/megahit/bin/megahit"  # 替换为 which megahit 的输出
export PATH="/home/gao/miniconda3/bin:$PATH"     # 添加环境变量

# 关键配置2：使用SSD临时目录
TMP_BASE="/ssd/tmp/megahit"
mkdir -p $TMP_BASE
chmod 1777 $TMP_BASE

# 生成样本列表
mapfile -t R1_FILES < <(find "$input_base" -name "*_paired_1.fastq.gz" | sort)
mapfile -t R2_FILES < <(find "$input_base" -name "*_paired_2.fastq.gz" | sort)
total=${#R1_FILES[@]}

# 并行控制
fifo="/tmp/megahit.fifo"
mkfifo $fifo || exit 1
exec 6<>$fifo
rm -f $fifo

for ((i=1; i<=max_parallel; i++)); do
    echo >&6
done

process_batch() {
    local batch_num=$1
    local start=$(( (batch_num-1)*batch_size ))
    local end=$(( batch_num*batch_size - 1 ))
    
    # 安全生成文件列表
    tmp_dir="$TMP_BASE/batch_$batch_num"
    mkdir -p $tmp_dir || { echo "无法创建临时目录"; exit 1; }
    
    # 关键修复：使用数组存储参数
    local r1_args=()
    local r2_args=()
    for ((i=start; i<=end && i<total; i++)); do
        r1_args+=("${R1_FILES[i]}")
        r2_args+=("${R2_FILES[i]}")
    done
    
    # 使用逗号分隔的字符串
    r1_list=$(IFS=,; echo "${r1_args[*]}")
    r2_list=$(IFS=,; echo "${r2_args[*]}")

    echo "启动批次 $batch_num (${#r1_args[@]} 个样本)"
    
    # 关键执行命令
    $MEGAHIT_PATH \
        -1 "$r1_list" \
        -2 "$r2_list" \
        -t $((threads/max_parallel)) \
        --continue \
        --memory 0.4 \
        --k-min 21 \
        --k-max 127 \
        --min-contig-len 500 \
        -o "${output_dir}/batch_${batch_num}" \
        --tmp-dir "$tmp_dir" || return $?

    # 清理
    rm -rf "$tmp_dir"
}

# 主流程
for ((batch=1; batch<=10; batch++)); do
    read -u6
    {
        process_batch $batch
        echo >&6
    } &
done

wait
exec 6>&-