#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# Config (EDIT if needed)
############################################
THREADS="${THREADS:-128}"

# repMAG 目录（dRep 输出）
REP_DIR="/home/gao/data/host_microbe/MAGs/drep/dereplicated_genomes"

# 输出目录
OUTDIR="/home/gao/data/host_microbe/MAGs/drep/checkm2_rerun"
mkdir -p "$OUTDIR"/{00_links,01_size,02_checkm2,log}

# CheckM2 docker + 数据库（按你现在的）
CHECKM2_IMAGE="quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1"
CHECKM2_DB_DMND="/home/gao/databases/CheckM2_database/uniref100.KO.1.dmnd"

############################################
# Checks
############################################
[[ -d "$REP_DIR" ]] || { echo "ERROR: REP_DIR not found: $REP_DIR" >&2; exit 1; }
[[ -s "$CHECKM2_DB_DMND" ]] || { echo "ERROR: CHECKM2_DB_DMND not found: $CHECKM2_DB_DMND" >&2; exit 1; }

############################################
# 0) Collect & link MAG FASTAs -> standardized .fna
############################################
echo "[0] Collect MAG fasta files from: $REP_DIR"

mapfile -t FASTAS < <(find "$REP_DIR" -maxdepth 1 -type f \( -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" \) | sort)
N=${#FASTAS[@]}
echo "  Found: $N"
[[ "$N" -gt 0 ]] || { echo "ERROR: no fasta found in $REP_DIR" >&2; exit 1; }

echo "[0.1] Stage MAG fastas as REAL files (*.fna) in: $OUTDIR/00_links"
rm -f "$OUTDIR/00_links/"*.fna 2>/dev/null || true

for f in "${FASTAS[@]}"; do
  bn="$(basename "$f")"
  mag="${bn%.*}"   # 去掉最后一个扩展名

  # 复制成真实文件，避免 docker 内软链接失效
  # rsync 比 cp 更稳（重复跑不会反复复制大文件）
  rsync -a --inplace "$f" "$OUTDIR/00_links/${mag}.fna"
done


LINK_N=$(ls -1 "$OUTDIR/00_links/"*.fna | wc -l | awk '{print $1}')
echo "  Linked: $LINK_N"

############################################
# 1) Compute genome size (bp / Mb)
############################################
echo "[1] Compute genome size (bp/Mb) ..."
SIZE_TSV="$OUTDIR/01_size/MAG_genome_size.tsv"

# awk 统计 fasta 总碱基数（去掉 header）
# 输出：MAG    genome_size_bp    genome_size_Mb
awk '
BEGIN{FS="\t"; OFS="\t"}
FNR==1{
  # file change hook not reliable in awk portable mode, so do in shell loop below
}
' /dev/null > /dev/null

: > "$SIZE_TSV"
echo -e "MAG\tgenome_size_bp\tgenome_size_Mb" > "$SIZE_TSV"

while IFS= read -r fa; do
  mag="$(basename "$fa")"
  mag="${mag%.fna}"

  # sum length of all sequence lines
  bp=$(awk 'BEGIN{n=0} !/^>/{gsub(/[ \r\t]/,""); n+=length($0)} END{print n}' "$fa")
  mb=$(awk -v bp="$bp" 'BEGIN{printf "%.6f", bp/1000000}')
  echo -e "${mag}\t${bp}\t${mb}" >> "$SIZE_TSV"
done < <(ls -1 "$OUTDIR/00_links/"*.fna | sort)

echo "  Wrote: $SIZE_TSV"

############################################
# 2) Run CheckM2 (docker)
############################################
echo "[2] Run CheckM2 (docker) ..."
CHECKM2_OUT="$OUTDIR/02_checkm2"

docker run --rm \
  -v "$OUTDIR/00_links":/input:ro \
  -v "$CHECKM2_OUT":/out \
  -v "$(dirname "$CHECKM2_DB_DMND")":"$(dirname "$CHECKM2_DB_DMND")":ro \
  "$CHECKM2_IMAGE" \
  bash -lc "
    set -euo pipefail
    checkm2 predict \
      --threads $THREADS \
      --input /input \
      --output-directory /out \
      --force \
      --database_path '$CHECKM2_DB_DMND' \
      -x fna
  " | tee "$OUTDIR/log/checkm2.predict.log"

# CheckM2 输出里通常有 quality_report.tsv
REPORT="$CHECKM2_OUT/quality_report.tsv"
if [[ ! -s "$REPORT" ]]; then
  echo "ERROR: CheckM2 report not found: $REPORT" >&2
  echo "Check log: $OUTDIR/log/checkm2.predict.log" >&2
  exit 1
fi
echo "  CheckM2 report: $REPORT"

############################################
# 3) Merge genome size + CheckM2 report
############################################
echo "[3] Merge genome size + CheckM2 report ..."
MERGED="$OUTDIR/MAG_genomeSize_checkm2.tsv"

python3 - <<'PY'
import pandas as pd, os

size_tsv = os.environ["SIZE_TSV"]
report   = os.environ["REPORT"]
out      = os.environ["MERGED"]

sz = pd.read_csv(size_tsv, sep="\t")
rp = pd.read_csv(report, sep="\t")

# CheckM2 第一列一般是 "name"（MAG名），但不同版本可能叫 genome/bin
id_candidates = [c for c in rp.columns if c.lower() in ["name","genome","bin_id","bin","mag","assembly"]]
idcol = id_candidates[0] if id_candidates else rp.columns[0]
rp = rp.rename(columns={idcol:"MAG"})

m = sz.merge(rp, on="MAG", how="left")
m.to_csv(out, sep="\t", index=False)
print("Wrote:", out, "rows=", m.shape[0], "cols=", m.shape[1])
PY \
  SIZE_TSV="$SIZE_TSV" \
  REPORT="$REPORT" \
  MERGED="$MERGED"

echo "[DONE]"
echo "Links:    $OUTDIR/00_links/*.fna"
echo "Size:     $SIZE_TSV"
echo "CheckM2:  $REPORT"
echo "Merged:   $MERGED"
