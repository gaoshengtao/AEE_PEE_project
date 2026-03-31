#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# rMAG proteins vs NCBI NR: protein similarity distribution
# 1) Prodigal predict proteins for each genome (.fna)
# 2) Concatenate proteins with MAG id in header
# 3) DIAMOND blastp to nr.dmnd (top1 hit)
# 4) Join with taxonomy (TSV/TXT) -> protein-level table
# 5) Summarize per class -> summary tables
# Output: tables only (no plots)
# ============================================================

# ---------- config ----------
GENOME_DIR="/home/gao/data/host_microbe/MAGs/drep/tree_input_genomes"
TAX_TSV="/home/gao/data/host_microbe/MAGs/Supplemental file 3-master_repMAG_taxonomy_by_rank.txt"
NR_DMND="/home/gao/backup/database/nr/nr.dmnd"
THREADS=64

OUTDIR="./nr_similarity_out"
mkdir -p "$OUTDIR"/{01_prodigal,02_diamond,03_tables,logs,tmp}

# ---------- checks ----------
command -v prodigal >/dev/null 2>&1 || { echo "[ERR] prodigal not found in PATH"; exit 1; }
command -v diamond  >/dev/null 2>&1 || { echo "[ERR] diamond not found in PATH"; exit 1; }
command -v python   >/dev/null 2>&1 || { echo "[ERR] python not found in PATH"; exit 1; }

[[ -d "$GENOME_DIR" ]] || { echo "[ERR] GENOME_DIR not found: $GENOME_DIR"; exit 1; }
[[ -s "$TAX_TSV" ]]    || { echo "[ERR] TAX_TSV not found: $TAX_TSV"; exit 1; }
[[ -s "$NR_DMND" ]]    || { echo "[ERR] NR_DMND not found: $NR_DMND"; exit 1; }

echo "[INFO] GENOME_DIR=$GENOME_DIR"
echo "[INFO] TAX_TSV=$TAX_TSV"
echo "[INFO] NR_DMND=$NR_DMND"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] THREADS=$THREADS"

# ---------- step1 outputs ----------
PROTEIN_DIR="$OUTDIR/01_prodigal"
ALLFAA="$OUTDIR/01_prodigal/all_rmag_proteins.faa"

# ---------- step2 outputs ----------
DIAMOND_OUT="$OUTDIR/02_diamond/rMAGs_vs_NR.top1.tsv"

# ---------- step3 outputs ----------
PROTEIN_TSV="$OUTDIR/03_tables/protein_similarity_with_taxonomy.tsv"
CLASS_SUMMARY="$OUTDIR/03_tables/summary_by_class.tsv"
CLASS_SUMMARY_ALL="$OUTDIR/03_tables/summary_by_class.ALL.tsv"

# ---------- step1: prodigal per genome ----------
# if ALLFAA exists and non-empty, skip step1 entirely
if [[ -s "$ALLFAA" ]]; then
  echo "[INFO] Step1 skip: found existing combined proteins: $ALLFAA"
else
  echo "[INFO] Step1: Prodigal predict proteins..."
  : > "$ALLFAA"
  n_genomes=0

  shopt -s nullglob
  for fna in "$GENOME_DIR"/*.fna; do
    [[ -s "$fna" ]] || continue
    n_genomes=$((n_genomes+1))
    base="$(basename "$fna" .fna)"
    faa="$PROTEIN_DIR/${base}.faa"

    if [[ -s "$faa" ]]; then
      echo "[INFO]   skip prodigal: $base (exists)"
    else
      echo "[INFO]   prodigal: $base"
      tmpfaa="$OUTDIR/tmp/${base}.raw.faa"
      prodigal -i "$fna" -a "$tmpfaa" -p meta -q 2>>"$OUTDIR/logs/prodigal.err.log"

      # prefix each protein header with MAG id: >MAGID|original_header
      awk -v mag="$base" 'BEGIN{OFS=""}
        /^>/ {sub(/^>/, ">"mag"|", $0)}
        {print}
      ' "$tmpfaa" > "$faa"
      rm -f "$tmpfaa"
    fi

    cat "$faa" >> "$ALLFAA"
  done
  shopt -u nullglob

  echo "[INFO] Step1 done. Genomes processed: $n_genomes"
  echo "[INFO] Combined proteins: $ALLFAA"
fi

# ---------- step2: diamond blastp top1 ----------
if [[ -s "$DIAMOND_OUT" ]]; then
  echo "[INFO] Step2 skip: found existing diamond result: $DIAMOND_OUT"
else
  echo "[INFO] Step2: DIAMOND blastp to NR (top1)..."
  diamond blastp \
    -d "$NR_DMND" \
    -q "$ALLFAA" \
    -o "$DIAMOND_OUT" \
    -f 6 qseqid sseqid pident length evalue bitscore qlen slen qcovhsp \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --query-cover 50 \
    --threads "$THREADS" \
    2> "$OUTDIR/logs/diamond.err.log"
  echo "[INFO] Step2 done: $DIAMOND_OUT"
fi

# ---------- step3: join with taxonomy (TSV) + summarize ----------
echo "[INFO] Step3: Join DIAMOND result with taxonomy (TSV) and summarize..."

python - <<PY
import csv, math
from collections import defaultdict

tax_tsv = r"""$TAX_TSV"""
diamond_tsv = r"""$DIAMOND_OUT"""
out_protein = r"""$PROTEIN_TSV"""
out_sum4 = r"""$CLASS_SUMMARY"""
out_sumall = r"""$CLASS_SUMMARY_ALL"""

# ---------- read taxonomy ----------
mag2class = {}
with open(tax_tsv, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    if reader.fieldnames is None:
        raise SystemExit("[ERR] TAX_TSV has no header row.")

    fields = set(reader.fieldnames)
    if "MAG_formal" not in fields:
        raise SystemExit(f"[ERR] TAX_TSV missing column: MAG_formal. Found: {reader.fieldnames}")

    if "class_final" in fields:
        class_col = "class_final"
    elif "Class" in fields:
        class_col = "Class"
    else:
        raise SystemExit(f"[ERR] TAX_TSV missing class column: need class_final or Class. Found: {reader.fieldnames}")

    for row in reader:
        mag = (row.get("MAG_formal") or "").strip()
        cls = (row.get(class_col) or "").strip() or "Others"
        if mag:
            mag2class[mag] = cls

if not mag2class:
    raise SystemExit("[ERR] No MAG_formal->class parsed from TAX_TSV (is the file tab-separated?)")

stats = defaultdict(lambda: {
    "n": 0,
    "sum_pident": 0.0,
    "sum_qcov": 0.0,
    "lt30": 0,
    "lt50": 0,
    "pidents": []
})

focus = {"Methanobacteria","Bacilli","Bacteroidia","Clostridia"}

def quantile(sorted_list, q):
    if not sorted_list:
        return float("nan")
    n = len(sorted_list)
    pos = (n - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(sorted_list[lo])
    return sorted_list[lo] + (sorted_list[hi] - sorted_list[lo]) * (pos - lo)

total_rows = 0
matched_rows = 0
kept_rows = 0

with open(diamond_tsv, "r", newline="") as fin, open(out_protein, "w", newline="") as fout:
    r = csv.reader(fin, delimiter="\\t")
    w = csv.writer(fout, delimiter="\\t")
    w.writerow(["qseqid","sseqid","MAG_formal","class_final","pident","length","evalue","bitscore","qlen","slen","qcovhsp"])

    for row in r:
        if len(row) < 9:
            continue
        total_rows += 1
        qseqid, sseqid = row[0], row[1]
        try:
            pident = float(row[2])
            alen   = float(row[3])
            evalue = float(row[4])
            bits   = float(row[5])
            qlen   = float(row[6])
            slen   = float(row[7])
            qcov   = float(row[8])
        except Exception:
            continue

        mag = qseqid.split("|", 1)[0]
        cls = mag2class.get(mag)
        if cls is None:
            continue
        matched_rows += 1

        if qcov < 50 or alen < 50:
            continue
        kept_rows += 1

        w.writerow([qseqid, sseqid, mag, cls, pident, int(alen), evalue, bits, int(qlen), int(slen), qcov])

        st = stats[cls]
        st["n"] += 1
        st["sum_pident"] += pident
        st["sum_qcov"] += qcov
        if pident < 30: st["lt30"] += 1
        if pident < 50: st["lt50"] += 1
        st["pidents"].append(pident)

def write_summary(path, only_focus=False):
    rows = []
    for cls, st in stats.items():
        if only_focus and cls not in focus:
            continue
        n = st["n"]
        p = sorted(st["pidents"])
        mean_p = (st["sum_pident"]/n) if n else float("nan")
        mean_q = (st["sum_qcov"]/n) if n else float("nan")
        rows.append({
            "class_final": cls,
            "n_proteins": n,
            "mean_pident": mean_p,
            "median_pident": quantile(p, 0.5),
            "q25_pident": quantile(p, 0.25),
            "q75_pident": quantile(p, 0.75),
            "frac_pident_lt30": (st["lt30"]/n) if n else float("nan"),
            "frac_pident_lt50": (st["lt50"]/n) if n else float("nan"),
            "mean_qcovhsp": mean_q
        })

    rows.sort(key=lambda x: (-x["n_proteins"], x["median_pident"] if not math.isnan(x["median_pident"]) else 1e9))

    with open(path, "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\\t")
            w.writeheader()
            w.writerows(rows)
        else:
            w = csv.writer(f, delimiter="\\t")
            w.writerow(["class_final","n_proteins","mean_pident","median_pident","q25_pident","q75_pident",
                        "frac_pident_lt30","frac_pident_lt50","mean_qcovhsp"])

write_summary(out_sumall, only_focus=False)
write_summary(out_sum4, only_focus=True)

print(f"[INFO] diamond rows total: {total_rows}")
print(f"[INFO] diamond rows matched to taxonomy: {matched_rows}")
print(f"[INFO] diamond rows kept after filters: {kept_rows}")
print("[INFO] protein-level table:", out_protein)
print("[INFO] class summary (ALL):", out_sumall)
print("[INFO] class summary (4 classes):", out_sum4)
PY

echo "[INFO] Step3 done."
echo "[INFO] Outputs:"
echo "  - $PROTEIN_TSV"
echo "  - $CLASS_SUMMARY_ALL"
echo "  - $CLASS_SUMMARY"

