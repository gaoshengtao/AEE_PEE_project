#!/usr/bin/env bash
# 利用 dRep Cdb + 全量 GTDB 注释，为 dRep 代表性 MAGs 继承 taxonomy，
# 并按 Bacteria / Archaea 分开构树（GTDB-Tk v2）

set -euo pipefail
shopt -s nullglob

########################################
# 0. 路径设置（按需修改）
########################################

# dRep 聚类结果
CDB="/home/gao/data/host_microbe/MAGs/drep/data_tables/Cdb.csv"

# 全部 MAG 的 GTDB 注释（之前 batch_run 跑过的结果）
BAC_SUM="/home/gao/data/host_microbe/MAGs/gtdbtk_out/batch_run/gtdbtk.bac120.summary.tsv"
ARC_SUM="/home/gao/data/host_microbe/MAGs/gtdbtk_out/batch_run/gtdbtk.ar53.summary.tsv"

# dRep 去冗余后的代表性 MAG 目录（已经统一为 .fna）
REP_DIR="/home/gao/data/host_microbe/MAGs/drep/tree_input_genomes"

# 代表 MAG 继承后的 taxonomy 输出表
OUT_TAX="/home/gao/data/host_microbe/MAGs/drep/drep_representatives.taxonomy.tsv"

# 按 domain 分开的代表性 MAG 目录
BAC_DIR="/home/gao/data/host_microbe/MAGs/drep/tree_input_bac"
ARC_DIR="/home/gao/data/host_microbe/MAGs/drep/tree_input_arc"

# 系统发育树输出目录
OUT_BAC="/home/gao/data/host_microbe/MAGs/drep/gtdbtk_tree_bac"
OUT_ARC="/home/gao/data/host_microbe/MAGs/drep/gtdbtk_tree_arc"

# 线程数
THREADS=32

echo "📌 dRep Cdb:         $CDB"
echo "📌 GTDB summary(bac): $BAC_SUM"
echo "📌 GTDB summary(arc): $ARC_SUM"
echo "📌 dRep 代表 MAG 目录: $REP_DIR"
echo "📌 继承 taxonomy 输出: $OUT_TAX"
echo

########################################
# 1. 检查 GTDB-Tk 数据库环境变量
########################################

if [[ -z "${GTDBTK_DATA_PATH:-}" ]]; then
    if [[ -d "$HOME/databases/gtdb" ]]; then
        export GTDBTK_DATA_PATH="$HOME/databases/gtdb"
        echo "ℹ️ 未检测到 GTDBTK_DATA_PATH，已自动设置为：$GTDBTK_DATA_PATH"
    else
        echo "❌ GTDBTK_DATA_PATH 未设置，且默认目录 ~/databases/gtdb 不存在，请手动配置。" >&2
        exit 1
    fi
fi

if ! command -v gtdbtk &>/dev/null; then
    echo "❌ 未找到 gtdbtk 命令，请在含 GTDB-Tk 的 conda 环境中运行此脚本。" >&2
    exit 1
fi

########################################
# 2. 用 Python 处理：
#    - 读 Cdb.csv + 两个 summary.tsv
#    - 为每个 cluster 选“最靠谱”的 taxonomy
#    - 为每个代表性 MAG 继承 taxonomy
#    - 按 Bacteria / Archaea 复制 fasta 到 BAC_DIR / ARC_DIR
########################################

export CDB BAC_SUM ARC_SUM REP_DIR OUT_TAX BAC_DIR ARC_DIR

python3 - << 'PY'
import csv
import os
import sys
import shutil

CDB      = os.environ["CDB"]
BAC_SUM  = os.environ["BAC_SUM"]
ARC_SUM  = os.environ["ARC_SUM"]
REP_DIR  = os.environ["REP_DIR"]
OUT_TAX  = os.environ["OUT_TAX"]
BAC_DIR  = os.environ["BAC_DIR"]
ARC_DIR  = os.environ["ARC_DIR"]

os.makedirs(BAC_DIR, exist_ok=True)
os.makedirs(ARC_DIR, exist_ok=True)

# 清空旧内容
for d in (BAC_DIR, ARC_DIR):
    for f in os.listdir(d):
        fp = os.path.join(d, f)
        if os.path.isfile(fp):
            os.remove(fp)

def load_taxonomy(summary_path):
    """
    读取 GTDB summary.tsv，返回:
    root_id -> taxonomy 字符串
    """
    tax = {}
    if not os.path.isfile(summary_path):
        return tax

    with open(summary_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return tax

        # 一般第一列是 user_genome，第二列是 gtdb_taxonomy
        for row in reader:
            if not row:
                continue
            genome_id = row[0].strip()
            if not genome_id:
                continue
            # 有的 summary 没给 taxonomy，就跳过
            if len(row) < 2:
                continue
            taxonomy = row[1].strip()
            root = genome_id  # 这里 summary 里一般是去掉后缀的名字，例如 14348 或 bin.1030
            tax[root] = taxonomy
    return tax

print("🔍 读取 GTDB 注释信息...", file=sys.stderr)
bac_tax = load_taxonomy(BAC_SUM)
arc_tax = load_taxonomy(ARC_SUM)

tax_by_root = {}
tax_by_root.update(bac_tax)
for k, v in arc_tax.items():
    # 如果同一个 root 同时出现在 bac & arc（基本不可能），优先保留以前的
    tax_by_root.setdefault(k, v)

print(f"   - 细菌注释条数: {len(bac_tax)}", file=sys.stderr)
print(f"   - 古菌注释条数: {len(arc_tax)}", file=sys.stderr)

# 读取 Cdb，建立：
#   root -> cluster_id
#   cluster_id -> [root1, root2, ...]
cluster_by_root = {}
roots_by_cluster = {}

print("🔍 读取 dRep Cdb.csv...", file=sys.stderr)
with open(CDB, "r", newline="") as f:
    reader = csv.DictReader(f)
    # 兼容不同列名：cluster 或 primary_cluster
    if "cluster" in reader.fieldnames:
        cluster_col = "cluster"
    elif "primary_cluster" in reader.fieldnames:
        cluster_col = "primary_cluster"
    else:
        raise SystemExit("Cdb.csv 中找不到 'cluster' 或 'primary_cluster' 列，请确认文件格式。")

    if "genome" not in reader.fieldnames:
        raise SystemExit("Cdb.csv 中找不到 'genome' 列，请确认文件格式。")

    for row in reader:
        genome = row["genome"].strip()
        if not genome:
            continue
        # basename 去掉路径 & 后缀
        base = os.path.basename(genome)
        root = base
        for suf in (".fa", ".fna", ".fasta"):
            if root.endswith(suf):
                root = root[:-len(suf)]
                break

        cluster_id = row[cluster_col].strip()
        if not cluster_id:
            continue
        cluster_by_root[root] = cluster_id
        roots_by_cluster.setdefault(cluster_id, []).append(root)

print(f"   - 共记录 cluster 数量: {len(roots_by_cluster)}", file=sys.stderr)

# 列出所有代表性 MAG（tree_input_genomes 下的 .fna）
rep_files = [f for f in os.listdir(REP_DIR) if f.endswith(".fna")]
root_to_repfile = {}
for f in rep_files:
    root = f
    for suf in (".fa", ".fna", ".fasta"):
        if root.endswith(suf):
            root = root[:-len(suf)]
            break
    root_to_repfile[root] = f

print(f"   - 代表性 MAG 文件数 (REP_DIR): {len(rep_files)}", file=sys.stderr)

# 为每个 cluster 选一个“最靠谱”的 taxonomy
def taxonomy_score(tx):
    """
    越长、unclassified 越少，得分越高
    """
    if not tx or tx.lower() == "na":
        return 0
    ranks = [r.strip() for r in tx.split(";") if r.strip()]
    score = 0
    for r in ranks:
        if "unclassified" in r or "norank" in r:
            continue
        score += 1
    return score

best_tax_by_cluster = {}
domain_by_cluster = {}

for cluster_id, roots in roots_by_cluster.items():
    best_tx = None
    best_score = -1
    best_domain = None

    for root in roots:
        tx = tax_by_root.get(root)
        if not tx:
            continue
        score = taxonomy_score(tx)
        if score > best_score:
            best_score = score
            best_tx = tx
            # 从 taxonomy 第一个 rank 提取 domain
            first = tx.split(";", 1)[0] if ";" in tx else tx
            first = first.strip()
            if first.startswith("d__Bacteria"):
                best_domain = "Bacteria"
            elif first.startswith("d__Archaea"):
                best_domain = "Archaea"
            else:
                best_domain = "Unknown"

    if best_tx:
        best_tax_by_cluster[cluster_id] = best_tx
        domain_by_cluster[cluster_id] = best_domain
    else:
        # 整个 cluster 都没任何 GTDB 注释
        best_tax_by_cluster[cluster_id] = "NA"
        domain_by_cluster[cluster_id] = "Unknown"

print("🔧 为每个 cluster 选定最佳 taxonomy 完成。", file=sys.stderr)

# 为每个代表性 root 继承 cluster 的 taxonomy，并复制到按 domain 的目录
rep_records = []

missing_cluster_roots = 0
missing_rep_file = 0

for root, rep_fname in root_to_repfile.items():
    cluster_id = cluster_by_root.get(root)
    if not cluster_id:
        missing_cluster_roots += 1
        print(f"⚠️ 代表 MAG {root} 不在 Cdb.csv 中，无法确定 cluster，跳过继承 taxonomy。", file=sys.stderr)
        continue

    tx = best_tax_by_cluster.get(cluster_id, "NA")
    dom = domain_by_cluster.get(cluster_id, "Unknown")

    # 代表 fasta 文件
    rep_path = os.path.join(REP_DIR, rep_fname)
    if not os.path.isfile(rep_path):
        missing_rep_file += 1
        print(f"⚠️ 代表 MAG 文件不存在: {rep_path}", file=sys.stderr)
        continue

    # 按 domain 复制 fasta
    if dom == "Bacteria":
        shutil.copy2(rep_path, os.path.join(BAC_DIR, rep_fname))
    elif dom == "Archaea":
        shutil.copy2(rep_path, os.path.join(ARC_DIR, rep_fname))
    else:
        # Unknown domain：先不复制，用 taxonomy 表里标出来
        print(f"⚠️ MAG {rep_fname} 继承的 domain=Unknown，暂不放入细菌/古菌树。", file=sys.stderr)

    rep_records.append((rep_fname, root, cluster_id, dom, tx))

print(f"📊 代表 MAG 中：", file=sys.stderr)
print(f"   - 未在 Cdb.csv 中找到 cluster 的: {missing_cluster_roots}", file=sys.stderr)
print(f"   - 找不到 fasta 文件的:           {missing_rep_file}", file=sys.stderr)

# 输出代表性 MAG 的 taxonomy 表
with open(OUT_TAX, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["genome_file", "root_id", "cluster_id", "domain", "gtdb_taxonomy_inherited"])
    for rec in rep_records:
        writer.writerow(rec)

print(f"✅ 已生成代表 MAG taxonomy 表: {OUT_TAX}", file=sys.stderr)
print(f"   - 细菌代表 fasta 数: {len(os.listdir(BAC_DIR))}", file=sys.stderr)
print(f"   - 古菌代表 fasta 数: {len(os.listdir(ARC_DIR))}", file=sys.stderr)

PY

echo
echo "✅ Python 阶段完成：已为代表 MAG 继承 taxonomy，并按 domain 分到了："
echo "   - Bacteria: $BAC_DIR"
echo "   - Archaea:  $ARC_DIR"
echo "   - 继承 taxonomy 表: $OUT_TAX"
echo

########################################
# 3. 用 GTDB-Tk 构建细菌树 / 古菌树
########################################

# 细菌树
if compgen -G "$BAC_DIR/*.fna" > /dev/null; then
    echo "🌳 构建细菌系统发育树 ..."
    mkdir -p "$OUT_BAC"

    gtdbtk de_novo_wf \
        --genome_dir "$BAC_DIR" \
        --out_dir "$OUT_BAC" \
        --bacteria \
        --outgroup_taxon p__Patescibacteriota \
        -x fna \
        --cpus "$THREADS"

    echo "✅ 细菌树完成，结果目录：$OUT_BAC"
else
    echo "⚠️ Bacteria 目录为空（$BAC_DIR），跳过细菌树构建。"
fi

echo

# 古菌树
if compgen -G "$ARC_DIR/*.fna" > /dev/null; then
    echo "🌲 构建古菌系统发育树 ..."
    mkdir -p "$OUT_ARC"

    gtdbtk de_novo_wf \
        --genome_dir "$ARC_DIR" \
        --out_dir "$OUT_ARC" \
        --archaea \
        --outgroup_taxon p__Altiarchaeota \
        -x fna \
        --cpus "$THREADS"

    echo "✅ 古菌树完成，结果目录：$OUT_ARC"
else
    echo "⚠️ Archaea 目录为空（$ARC_DIR），跳过古菌树构建。"
fi

echo
echo "🎉 全流程结束："
echo "   - 代表 MAG taxonomy: $OUT_TAX"
echo "   - 细菌树:            $OUT_BAC"
echo "   - 古菌树:            $OUT_ARC"

