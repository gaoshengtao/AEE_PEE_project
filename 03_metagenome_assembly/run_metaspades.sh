#!/bin/bash
# 宏基因组组装并行处理脚本 v3.2（带完成检测）

# >>> 配置区 <<<
input_root="/home/gao/data/host_microbe/hostfree/"
output_root="/home/gao/data/host_microbe/assembly/metaspades/"
nanopore_file="/home/gao/data/host_microbe/hostfree/ont.filtered.combined.fastq.gz"
spades_exec="/home/gao/apps/SPAdes-4.0.0/bin/metaspades.py"
threads=96
memory=256
max_jobs=6
log_dir="${output_root}/logs"
interval_check=120

# >>> 初始化 <<<
mkdir -p "$output_root" "$log_dir"
declare -i running=0 completed=0 skipped=0 failed=0
timestamp=$(date +%Y%m%d_%H%M)
echo "====== 批量处理开始 [$timestamp] ======" > "${log_dir}/master.log"

# 生成样本队列
mapfile -t sample_list < <(find "$input_root" -maxdepth 1 -type d -name "H*" -exec basename {} \; | sort)

# >>> 任务调度引擎 <<<
for sample_id in "${sample_list[@]}"; do
    output_dir="${output_root}/${sample_id}"
    done_marker1="${output_dir}/contigs.fasta"
    done_marker2="${output_dir}/scaffolds.fasta"
    
    # 跳过已完成的样本（新增双文件校验）
    if [[ -f "$done_marker1" && -f "$done_marker2" ]]; then
        echo "[跳过] $sample_id 已存在完整结果" | tee -a "${log_dir}/skipped.log"
        ((skipped++))
        continue
    fi

    # 输入验证（保持原逻辑）
    r1="${input_root}/${sample_id}/${sample_id}_paired_1.fastq.gz"
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "[错误] $sample_id 输入文件缺失" | tee -a "${log_dir}/error.log"
        ((failed++))
        continue
    fi

    # 创建任务目录
    mkdir -p "$output_dir"
    
    # 生成并行任务（保持原队列逻辑）
    (
        echo "[启动][$(date)] $sample_id 开始处理" | tee -a "${log_dir}/${sample_id}.log"
        "$spades_exec" \
            -1 "$r1" \
            -2 "$r2" \
            --nanopore "$nanopore_file" \
            -o "$output_dir" \
            -t "$threads" \
            -m "$memory" \
            >> "${log_dir}/${sample_id}.log" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[完成][$(date)] $sample_id" | tee -a "${log_dir}/success.log"
            ((completed++))
        else
            echo "[失败] $sample_id" | tee -a "${log_dir}/error.log"
            ((failed++))
        fi
    ) &
    
    # 控制并行度（保持原逻辑）
    if (( $(jobs -rp | wc -l) >= max_jobs )); then
        wait -n
    fi
done

# 等待剩余任务
wait

# >>> 生成报告 <<<
echo "====== 处理完成 [$(date)] ======" | tee -a "${log_dir}/master.log"
echo "成功: $completed | 跳过: $skipped | 失败: $failed" | tee -a "${log_dir}/master.log"