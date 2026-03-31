#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# User config (edit if needed)
############################################
DREP_DIR="/home/gao/data/host_microbe/MAGs/drep/dereplicated_genomes"
WORK_BASE="/home/gao/data/host_microbe/MAGs/Phylotree/rebuild_$(date +%Y%m%d_%H%M%S)"

# PhyloPhlAn database root (contains phylophlan.dmnd / phylophlan.faa)
PHYLO_DB_ROOT="/home/gao/data/host_microbe/MAGs/Phylotree/phylophlan_databases"
DB_NAME="phylophlan"

# Which config to use
CFG_DIR="/home/gao/data/host_microbe/MAGs/Phylotree"
CFG="supertree_aa.cfg"   # 你目录里有这个文件

# Threads
NPROC="${NPROC:-128}"

############################################
# Derived paths
############################################
IN_CLEAN="${WORK_BASE}/01_input_clean_fna"
OUT_DIR="${WORK_BASE}/02_phylophlan_out"
LOG_DIR="${WORK_BASE}/logs"
MANIFEST="${WORK_BASE}/input_manifest.tsv"

mkdir -p "$IN_CLEAN" "$OUT_DIR" "$LOG_DIR"

echo "[INFO] DREP_DIR     = $DREP_DIR"
echo "[INFO] WORK_BASE    = $WORK_BASE"
echo "[INFO] IN_CLEAN     = $IN_CLEAN"
echo "[INFO] OUT_DIR      = $OUT_DIR"
echo "[INFO] PHYLO_DB_ROOT= $PHYLO_DB_ROOT"
echo "[INFO] DB_NAME      = $DB_NAME"
echo "[INFO] CFG          = $CFG_DIR/$CFG"
echo "[INFO] NPROC        = $NPROC"

############################################
# 0) Basic checks
############################################
command -v phylophlan >/dev/null 2>&1 || { echo "[ERROR] phylophlan not in PATH"; exit 1; }

if [ ! -d "$DREP_DIR" ]; then
  echo "[ERROR] DREP_DIR not found: $DREP_DIR" >&2
  exit 1
fi

if [ ! -f "${CFG_DIR}/${CFG}" ]; then
  echo "[ERROR] config not found: ${CFG_DIR}/${CFG}" >&2
  exit 1
fi

if [ ! -f "${PHYLO_DB_ROOT}/${DB_NAME}/${DB_NAME}.dmnd" ] && [ ! -f "${PHYLO_DB_ROOT}/${DB_NAME}/${DB_NAME}.faa" ]; then
  # 兼容你 tree 显示的路径结构：phylophlan_databases/phylophlan/phylophlan.dmnd
  if [ -f "${PHYLO_DB_ROOT}/${DB_NAME}/${DB_NAME}.dmnd" ]; then
    :
  fi
fi

# Try to locate dmnd/faa robustly
DB_DMND="$(ls -1 ${PHYLO_DB_ROOT}/${DB_NAME}/${DB_NAME}.dmnd 2>/dev/null || true)"
DB_FAA="$(ls -1 ${PHYLO_DB_ROOT}/${DB_NAME}/${DB_NAME}.faa  2>/dev/null || true)"
if [ -z "$DB_DMND" ] && [ -z "$DB_FAA" ]; then
  # fallback: search within root
  DB_DMND="$(find "$PHYLO_DB_ROOT" -maxdepth 3 -type f -name "${DB_NAME}.dmnd" | head -n 1 || true)"
  DB_FAA="$(find "$PHYLO_DB_ROOT" -maxdepth 3 -type f -name "${DB_NAME}.faa"  | head -n 1 || true)"
fi
if [ -z "$DB_DMND" ] && [ -z "$DB_FAA" ]; then
  echo "[ERROR] Cannot find ${DB_NAME}.dmnd or ${DB_NAME}.faa under: $PHYLO_DB_ROOT" >&2
  exit 1
fi
echo "[INFO] DB_DMND = ${DB_DMND:-NA}"
echo "[INFO] DB_FAA  = ${DB_FAA:-NA}"

############################################
# 1) Prepare clean input genomes
#    - unify extension to .fna
#    - sanitize headers (remove spaces etc.)
#    - ensure unique IDs (detect duplicates)
############################################
echo "==[1] Prepare clean genomes (.fna) + safe headers =="

# collect input files
mapfile -t GENOMES < <(find "$DREP_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | sort)

if [ "${#GENOMES[@]}" -eq 0 ]; then
  echo "[ERROR] No genome files found in $DREP_DIR" >&2
  exit 1
fi

echo -e "orig_path\tclean_id\tclean_path\tbytes" > "$MANIFEST"

# helper: make safe id from filename (strip extension)
make_id() {
  local f="$1"
  local base
  base="$(basename "$f")"
  base="${base%.*}"          # strip last extension
  # keep only safe chars (letters/digits/._-)
  base="$(echo "$base" | tr ' ' '_' | sed 's/[^A-Za-z0-9._-]/_/g')"
  echo "$base"
}

# detect duplicate IDs
tmp_ids="${WORK_BASE}/tmp.ids"
: > "$tmp_ids"
for f in "${GENOMES[@]}"; do
  make_id "$f" >> "$tmp_ids"
done

dups="$(sort "$tmp_ids" | uniq -d | head -n 20 || true)"
if [ -n "$dups" ]; then
  echo "[ERROR] Duplicate genome IDs detected (showing up to 20):"
  echo "$dups"
  echo "Fix duplicates in filenames first (dRep output sometimes has _dup1 etc.)." >&2
  exit 1
fi

# copy + sanitize FASTA headers
for f in "${GENOMES[@]}"; do
  id="$(make_id "$f")"
  out="${IN_CLEAN}/${id}.fna"

  # sanitize:
  #  - keep sequence as-is
  #  - rewrite header to: >ID  (single token)
  # This avoids "spaces in header" breaking downstream sample naming.
  awk -v id="$id" '
    BEGIN{OFS=""}
    /^>/ {print ">", id; next}
    {gsub(/\r/,""); print}
  ' "$f" > "$out"

  bytes="$(wc -c < "$out" | tr -d ' ')"
  echo -e "${f}\t${id}\t${out}\t${bytes}" >> "$MANIFEST"
done

echo "[INFO] Input genomes found      : ${#GENOMES[@]}"
echo "[INFO] Clean genomes generated  : $(ls -1 "$IN_CLEAN"/*.fna | wc -l)"

############################################
# 2) Run PhyloPhlAn
############################################
echo "==[2] Run PhyloPhlAn (AA supertree) =="

# NOTE:
# - --genome_extension fna ensures it only picks .fna in IN_CLEAN
# - -d phylophlan uses DB_NAME
# - --databases_folder points to PHYLO_DB_ROOT
# - --configs_folder points to CFG_DIR, and --config_file is CFG
# - --nproc for speed
# - --diversity high is safer for mixed MAGs; you can change to low/medium if needed
phylophlan \
  -i "$IN_CLEAN" \
  -o "$OUT_DIR" \
  -d "$DB_NAME" \
  --databases_folder "$PHYLO_DB_ROOT" \
  --configs_folder "$CFG_DIR" \
  --config_file "$CFG" \
  --genome_extension "fna" \
  --nproc "$NPROC" \
  --diversity "high" \
  --verbose 2>&1 | tee "${LOG_DIR}/phylophlan.run.log"

############################################
# 3) Quick diagnostics: how many made it into concatenated alignment / tree tips
############################################
echo "==[3] Diagnostics =="

echo "[INFO] Expected (clean input) genomes: $(ls -1 "$IN_CLEAN"/*.fna | wc -l)"

# Find concatenated alignment (PhyloPhlAn output varies by config; search common patterns)
ALN="$(find "$OUT_DIR" -type f \( -name "*concatenated*.aln" -o -name "*supermatrix*.aln" -o -name "*.aln" \) | head -n 1 || true)"
if [ -n "$ALN" ]; then
  echo "[INFO] Alignment file: $ALN"
  nseq="$(grep -c '^>' "$ALN" || true)"
  echo "[INFO] Sequences in alignment: $nseq"
else
  echo "[WARN] Cannot find alignment (*.aln) in $OUT_DIR"
fi

TREE="$(find "$OUT_DIR" -type f \( -name "*.treefile" -o -name "*.tree" -o -name "*.nwk" \) | head -n 1 || true)"
if [ -n "$TREE" ]; then
  echo "[INFO] Tree file: $TREE"
  # count tips roughly (works for Newick without branch labels weirdness)
  # better: use R/ete3 later; here is quick approx:
  tips_est="$(python3 - <<'PY'
import re,sys
p=sys.argv[1]
s=open(p,'r',encoding='utf-8',errors='ignore').read().strip()
# rough: split by commas then count leaf labels before : or ) or ;
# Not perfect but good sanity check
leaves=re.findall(r'([^\(\),:;]+)\s*:', s)
print(len(leaves))
PY
"$TREE")"
  echo "[INFO] Estimated tips in tree : $tips_est"
else
  echo "[WARN] Cannot find tree file in $OUT_DIR"
fi

echo "[DONE] All outputs under: $WORK_BASE"

