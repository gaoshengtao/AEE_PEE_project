#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# Inputs you already have
############################################
GENOME_DIR="${1:-/home/gao/data/host_microbe/MAGs/drep/tree_input_genomes}"

# where to write everything
OUT_DIR="${2:-/home/gao/data/host_microbe/MAGs/drep/phylophlan_tree_from_prepared}"

THREADS="${THREADS:-64}"

# If you already have predicted proteins, point here to skip prodigal.
# Example: export PROTEIN_DIR=/path/to/existing_faa_dir
PROTEIN_DIR="${PROTEIN_DIR:-}"

############################################
# Real PhyloPhlAn DB paths (ABSOLUTE)
############################################
PHYLO_BASE="/home/gao/data/host_microbe/MAGs/Phylotree"
DB_REAL="$(readlink -f "$PHYLO_BASE/phylophlan_databases")"
DMND_REAL="$DB_REAL/phylophlan/phylophlan.dmnd"

############################################
# Checks
############################################
command -v phylophlan >/dev/null 2>&1 || { echo "ERROR: phylophlan not in PATH"; exit 1; }
command -v prodigal  >/dev/null 2>&1 || { echo "ERROR: prodigal not in PATH"; exit 1; }
[ -d "$GENOME_DIR" ] || { echo "ERROR: GENOME_DIR not found: $GENOME_DIR" >&2; exit 1; }
[ -f "$DMND_REAL" ] || { echo "ERROR: DMND not found: $DMND_REAL" >&2; exit 1; }

mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/work"
PROTS="$WORK/01_proteins_faa"
LOG="$OUT_DIR/run.log"
mkdir -p "$WORK" "$PROTS"

echo "========== [PhyloPhlAn | prepared genomes | ABS DB in cfg] ==========" | tee "$LOG"
echo "[INFO] GENOME_DIR = $GENOME_DIR" | tee -a "$LOG"
echo "[INFO] OUT_DIR    = $OUT_DIR" | tee -a "$LOG"
echo "[INFO] THREADS    = $THREADS" | tee -a "$LOG"
echo "[INFO] DB_REAL    = $DB_REAL" | tee -a "$LOG"
echo "[INFO] DMND_REAL  = $DMND_REAL" | tee -a "$LOG"
echo "[INFO] PROTEIN_DIR= ${PROTEIN_DIR:-<none>}" | tee -a "$LOG"

n_genomes=$(ls -1 "$GENOME_DIR"/*.fna 2>/dev/null | wc -l | awk '{print $1}')
echo "[INFO] N genomes (*.fna) = $n_genomes" | tee -a "$LOG"
[ "$n_genomes" -gt 0 ] || { echo "ERROR: no .fna in $GENOME_DIR" | tee -a "$LOG" >&2; exit 1; }

############################################
# 1) Prepare proteins (reuse if provided)
############################################
if [ -n "${PROTEIN_DIR:-}" ]; then
  [ -d "$PROTEIN_DIR" ] || { echo "ERROR: PROTEIN_DIR not found: $PROTEIN_DIR" | tee -a "$LOG" >&2; exit 1; }
  echo "[INFO] Using existing proteins from: $PROTEIN_DIR" | tee -a "$LOG"
  PROT_INPUT="$(readlink -f "$PROTEIN_DIR")"
else
  echo "[INFO] Predict proteins with prodigal -> $PROTS" | tee -a "$LOG"
  rm -f "$PROTS"/*.faa 2>/dev/null || true

  ls -1 "$GENOME_DIR"/*.fna | \
    xargs -I{} -P "$THREADS" bash -c '
      g="{}"
      id=$(basename "$g" .fna)
      out="'"$PROTS"'/${id}.faa"
      prodigal -i "$g" -a "$out" -p meta -q >/dev/null 2>&1 || {
        echo "[WARN] prodigal failed: $g" >&2
        rm -f "$out"
      }
    '

  n_faa=$(ls -1 "$PROTS"/*.faa 2>/dev/null | wc -l | awk '{print $1}')
  echo "[INFO] Proteins predicted = $n_faa" | tee -a "$LOG"
  [ "$n_faa" -gt 0 ] || { echo "ERROR: no faa produced" | tee -a "$LOG" >&2; exit 1; }
  PROT_INPUT="$(readlink -f "$PROTS")"
fi

############################################
# 3) Run phylophlan (AA supermatrix)
############################################
mkdir -p "$OUT_DIR/phylophlan_out"
cd "$OUT_DIR"

echo "[INFO] Running phylophlan..." | tee -a "$LOG"
phylophlan \
  -i "$PROT_INPUT" \
  -o "$OUT_DIR/phylophlan_out" \
  -d phylophlan \
  --data_folder "$DB_REAL" \
  --diversity high \
  --fast \
  --nproc "$THREADS" \
  --config_file /home/gao/data/host_microbe/MAGs/supermatrix_aa.cfg \
  2>&1 | tee -a "$LOG"

############################################
# 4) Collect outputs
############################################
echo "[INFO] Collecting outputs..." | tee -a "$LOG"

ALN=$(find "$OUT_DIR/phylophlan_out" -type f \( -name "*concatenated*.aln" -o -name "*.aln" \) | head -n 1 || true)
if [ -n "${ALN:-}" ] && [ -f "$ALN" ]; then
  echo "[INFO] Alignment file: $ALN" | tee -a "$LOG"
  IN_ALN_IDS="$OUT_DIR/in_alignment_ids.txt"
  grep '^>' "$ALN" | sed 's/^>//; s/[[:space:]].*$//' | sort -u > "$IN_ALN_IDS"
  echo "[INFO] N genomes in alignment = $(wc -l "$IN_ALN_IDS" | awk '{print $1}')" | tee -a "$LOG"
else
  echo "[WARN] No alignment found; see $LOG" | tee -a "$LOG"
fi

TREE_CAND="$OUT_DIR/tree_candidates.txt"
find "$OUT_DIR/phylophlan_out" -type f \( -name "*.tre" -o -name "*.tree" -o -name "*.nwk" -o -name "*bestTree*" -o -name "*final_tree*" \) \
  | sort -u > "$TREE_CAND" || true
if [ "$(wc -l "$TREE_CAND" | awk '{print $1}')" -gt 0 ]; then
  FINAL_TREE="$OUT_DIR/final_tree.nwk"
  cp -f "$(head -n 1 "$TREE_CAND")" "$FINAL_TREE"
  echo "[INFO] Final tree: $FINAL_TREE" | tee -a "$LOG"
else
  echo "[WARN] No tree detected automatically; see outputs." | tee -a "$LOG"
fi

echo "[DONE] See log: $LOG" | tee -a "$LOG"

