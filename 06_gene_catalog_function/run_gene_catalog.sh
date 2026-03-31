#!/usr/bin/env bash
set -euo pipefail

############################################
# Pipeline: ORF -> NR Gene Catalog -> eggNOG -> kallisto quant
############################################

# ---------------------------
# 0) User-configurable inputs
# ---------------------------
MAG_DIR="/home/gao/data/host_microbe/MAGs/drep/dereplicated_genomes"
HOSTFREE_DIR="/home/gao/data/host_microbe/hostfree"

# Output root
OUTDIR="/home/gao/data/host_microbe/gene_catalog"
mkdir -p "$OUTDIR"

# Threads
THREADS="${THREADS:-32}"

# cd-hit-est / kallisto are in metagenome conda env
METAGENOME_ENV="${METAGENOME_ENV:-metagenome}"

# cd-hit-est parameters (edit if needed)
CDHIT_ID="${CDHIT_ID:-0.95}"     # sequence identity
CDHIT_A_S="${CDHIT_A_S:-0.9}"    # alignment coverage for shorter sequence
CDHIT_M="${CDHIT_M:-0}"          # memory limit in MB (0 = unlimited)

# eggNOG-mapper (as you provided)
EGGNOG_MAPPER="$HOME/apps/eggnog-mapper-2.1.13/emapper.py"
EGGNOG_DB_DIR="/home/gao/databases/eggnog"
EGGNOG_CONDA_ENV="${EGGNOG_CONDA_ENV:-eggnog}"
EMAPPER_CPU="${EMAPPER_CPU:-8}"

# kallisto settings
KALLISTO_BOOTSTRAPS="${KALLISTO_BOOTSTRAPS:-0}"   # 0 means no bootstraps
KALLISTO_FRAGLEN="${KALLISTO_FRAGLEN:-200}"       # only used for single-end; keep for completeness
KALLISTO_SD="${KALLISTO_SD:-20}"                  # only used for single-end

# ---------------------------
# 1) Derived paths
# ---------------------------
LOGDIR="$OUTDIR/logs"
ORFDIR="$OUTDIR/01_orf"
CATDIR="$OUTDIR/02_catalog"
NRDIR="$OUTDIR/03_NR_cdhitest"
EGGDIR="$OUTDIR/04_eggnog"
KALDIR="$OUTDIR/05_kallisto"

mkdir -p "$LOGDIR" "$ORFDIR" "$CATDIR" "$NRDIR" "$EGGDIR" "$KALDIR"

ALL_GENES_FNA="$CATDIR/all_genes.fna"
ALL_PROTS_FAA="$CATDIR/all_proteins.faa"

NR_GENES_FNA="$NRDIR/NR_genes.fna"
NR_PROTS_FAA="$NRDIR/NR_proteins.faa"
NR_CDHIT_CLSTR="$NRDIR/NR_genes.fna.clstr"

KALLISTO_INDEX="$KALDIR/NR_genes.kallisto.idx"

# ---------------------------
# 2) Helpers
# ---------------------------
timestamp() { date +"%F %T"; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOGDIR/pipeline.log"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Missing command: $1" >&2; exit 1; }
}

conda_run() {
  # usage: conda_run <env> <cmd...>
  local env="$1"; shift
  conda run -n "$env" "$@"
}

# ---------------------------
# 3) Check dependencies
# ---------------------------
log "Checking dependencies..."
need_cmd awk
need_cmd sed
need_cmd grep
need_cmd find
need_cmd gzip
need_cmd xargs
need_cmd conda


# eggnog env check
conda_run "$EGGNOG_CONDA_ENV" python -V >/dev/null 2>&1 || {
  echo "[ERROR] eggnog conda env not usable: $EGGNOG_CONDA_ENV" >&2; exit 1;
}
[[ -x "$EGGNOG_MAPPER" ]] || { echo "[ERROR] EGGNOG_MAPPER not executable: $EGGNOG_MAPPER" >&2; exit 1; }
[[ -d "$EGGNOG_DB_DIR" ]] || { echo "[ERROR] EGGNOG_DB_DIR not found: $EGGNOG_DB_DIR" >&2; exit 1; }

log "All dependencies look OK."

# ---------------------------
# 4) Collect MAG fasta files
# ---------------------------
log "Scanning MAGs in: $MAG_DIR"
mapfile -t MAG_FASTAS < <(find "$MAG_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | sort)

if [[ "${#MAG_FASTAS[@]}" -eq 0 ]]; then
  echo "[ERROR] No MAG fasta found in: $MAG_DIR" >&2
  exit 1
fi
log "Found ${#MAG_FASTAS[@]} MAG genomes."

# ---------------------------
# 5) ORF prediction per MAG (Prokka not requested; use Prodigal)
#    Output:
#      <MAG>.genes.fna (nucleotides)
#      <MAG>.proteins.faa (proteins)
# ---------------------------
log "Step 1: ORF prediction with Prodigal (-p meta)..."
for f in "${MAG_FASTAS[@]}"; do
  base="$(basename "$f")"
  mag="${base%.*}"

  out_fna="$ORFDIR/${mag}.genes.fna"
  out_faa="$ORFDIR/${mag}.proteins.faa"
  out_gff="$ORFDIR/${mag}.prodigal.gff"
  out_sco="$ORFDIR/${mag}.prodigal.sco"

  if [[ -s "$out_fna" && -s "$out_faa" ]]; then
    log "  [skip] $mag ORFs already exist."
    continue
  fi

  log "  [run] $mag"
  conda_run "$METAGENOME_ENV" prodigal \
    -i "$f" \
    -p meta \
    -d "$out_fna" \
    -a "$out_faa" \
    -f gff \
    -o "$out_gff" \
    -s "$out_sco" \
    2>>"$LOGDIR/prodigal.stderr.log"

  # Prefix headers to keep uniqueness across MAGs
  # Prodigal headers like: >contig_1_1 # 1 # 123 # ...
  # We convert to: >MAG|contig_1_1_1  (keep only first token)
  tmp_fna="${out_fna}.tmp"
  tmp_faa="${out_faa}.tmp"
  awk -v PFX="${mag}|" 'BEGIN{FS=" "} /^>/{sub(/^>/,">"PFX,$1); print $1; next} {print}' "$out_fna" > "$tmp_fna"
  awk -v PFX="${mag}|" 'BEGIN{FS=" "} /^>/{sub(/^>/,">"PFX,$1); print $1; next} {print}' "$out_faa" > "$tmp_faa"
  mv "$tmp_fna" "$out_fna"
  mv "$tmp_faa" "$out_faa"
done
log "Step 1 done."

# ---------------------------
# 6) Merge all genes/proteins
# ---------------------------
log "Step 2: Merge all genes/proteins to build raw catalog..."
: > "$ALL_GENES_FNA"
: > "$ALL_PROTS_FAA"

# concatenate in sorted order for reproducibility
find "$ORFDIR" -type f -name "*.genes.fna" | sort | while read -r f; do
  cat "$f" >> "$ALL_GENES_FNA"
done
find "$ORFDIR" -type f -name "*.proteins.faa" | sort | while read -r f; do
  cat "$f" >> "$ALL_PROTS_FAA"
done

log "Raw catalog created:"
log "  - $ALL_GENES_FNA"
log "  - $ALL_PROTS_FAA"

# ---------------------------
# 7) cd-hit-est dereplication to build NR gene catalog
#    Input: all_genes.fna
#    Output: NR_genes.fna + cluster file
# ---------------------------
log "Step 3: cd-hit-est dereplication..."
if [[ -s "$NR_GENES_FNA" && -s "$NR_CDHIT_CLSTR" ]]; then
  log "  [skip] NR catalog already exists."
else
  conda_run "$METAGENOME_ENV" cd-hit-est \
    -i "$ALL_GENES_FNA" \
    -o "$NR_GENES_FNA" \
    -c "$CDHIT_ID" \
    -aS "$CDHIT_A_S" \
    -T "$THREADS" \
    -M "$CDHIT_M" \
    -d 0 \
    > "$LOGDIR/cdhitest.stdout.log" 2> "$LOGDIR/cdhitest.stderr.log"

  log "NR genes written:"
  log "  - $NR_GENES_FNA"
  log "  - $NR_CDHIT_CLSTR"
fi

# ---------------------------
# 8) Build NR proteins corresponding to NR genes
#    Approach:
#      - use cd-hit-est output representatives (NR_genes.fna) headers
#      - fetch matching proteins from ALL_PROTS_FAA by ID
#
#    We assume Prodigal produced gene/protein IDs consistent:
#      nucleotide header == protein header (after our prefixing)
# ---------------------------
log "Step 4: Extract NR proteins for NR genes..."

if [[ -s "$NR_PROTS_FAA" ]]; then
  log "  [skip] NR proteins already exist."
else
  # Build a set of representative IDs
  REP_IDS="$NRDIR/NR_gene_ids.txt"
  grep '^>' "$NR_GENES_FNA" | sed 's/^>//' > "$REP_IDS"

  # FASTA fetch without seqkit: use awk two-pass
  # First pass: load IDs into hash; second pass: print matched records
  awk 'NR==FNR{ids[$1]=1; next}
       /^>/{h=substr($0,2); keep=(h in ids)}
       {if(keep) print}
  ' "$REP_IDS" "$ALL_PROTS_FAA" > "$NR_PROTS_FAA"

  # sanity
  nr_n=$(grep -c '^>' "$NR_GENES_FNA" || true)
  nr_p=$(grep -c '^>' "$NR_PROTS_FAA" || true)
  log "  NR genes:     $nr_n"
  log "  NR proteins:  $nr_p"
  if [[ "$nr_n" -ne "$nr_p" ]]; then
    log "[WARN] NR gene count != NR protein count. This usually means header mismatch between .fna and .faa."
    log "       Check a few IDs from $NRDIR/NR_gene_ids.txt vs $ALL_PROTS_FAA."
  fi
fi

# ---------------------------
# 9) eggNOG functional annotation (on NR proteins)
# ---------------------------
log "Step 5: eggNOG-mapper annotation..."
EGG_PREFIX="$EGGDIR/NR"
EGG_DONE_FLAG="$EGGDIR/.eggnog.done"

if [[ -f "$EGG_DONE_FLAG" ]]; then
  log "  [skip] eggNOG step already done."
else
  mkdir -p "$EGGDIR"
  conda_run "$EGGNOG_CONDA_ENV" python "$EGGNOG_MAPPER" \
    -i "$NR_PROTS_FAA" \
    --output "$EGG_PREFIX" \
    --output_dir "$EGGDIR" \
    --data_dir "$EGGNOG_DB_DIR" \
    --cpu "$EMAPPER_CPU" \
    --itype proteins \
    --override \
    > "$LOGDIR/eggnog.stdout.log" 2> "$LOGDIR/eggnog.stderr.log"

  # Mark done
  touch "$EGG_DONE_FLAG"
  log "eggNOG outputs in: $EGGDIR"
  log "Typical key file: $EGGDIR/NR.emapper.annotations"
fi

# ---------------------------
# 10) kallisto quantification
#      - build index from NR genes (as transcripts)
#      - run quant per sample (paired fastq)
# ---------------------------
log "Step 6: kallisto index..."
if [[ -s "$KALLISTO_INDEX" ]]; then
  log "  [skip] kallisto index exists: $KALLISTO_INDEX"
else
  conda_run "$METAGENOME_ENV" kallisto index \
    -i "$KALLISTO_INDEX" \
    "$NR_GENES_FNA" \
    > "$LOGDIR/kallisto_index.stdout.log" 2> "$LOGDIR/kallisto_index.stderr.log"
  log "kallisto index built: $KALLISTO_INDEX"
fi

log "Step 7: kallisto quant for each sample under: $HOSTFREE_DIR"

# Discover paired fastqs: HOSTFREE_DIR/<Sample>/<Sample>_paired_1.fastq.gz
mapfile -t R1_LIST < <(find "$HOSTFREE_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_paired_1.fastq.gz" | sort)

if [[ "${#R1_LIST[@]}" -eq 0 ]]; then
  log "[ERROR] No R1 fastqs found like *_paired_1.fastq.gz in $HOSTFREE_DIR/*/"
  exit 1
fi

log "Found ${#R1_LIST[@]} samples with paired reads."

for r1 in "${R1_LIST[@]}"; do
  r2="${r1/_paired_1.fastq.gz/_paired_2.fastq.gz}"
  if [[ ! -s "$r2" ]]; then
    log "[WARN] Missing R2 for: $r1 ; skip."
    continue
  fi

  sample_dir="$(dirname "$r1")"
  sample="$(basename "$sample_dir")"   # e.g. H0619
  out="$KALDIR/$sample"
  mkdir -p "$out"

  if [[ -s "$out/abundance.tsv" ]]; then
    log "  [skip] $sample kallisto done."
    continue
  fi

  log "  [run] $sample"
  if [[ "$KALLISTO_BOOTSTRAPS" -gt 0 ]]; then
    conda_run "$METAGENOME_ENV" kallisto quant \
      -i "$KALLISTO_INDEX" \
      -o "$out" \
      -t "$THREADS" \
      -b "$KALLISTO_BOOTSTRAPS" \
      "$r1" "$r2" \
      > "$LOGDIR/kallisto_${sample}.stdout.log" 2> "$LOGDIR/kallisto_${sample}.stderr.log"
  else
    conda_run "$METAGENOME_ENV" kallisto quant \
      -i "$KALLISTO_INDEX" \
      -o "$out" \
      -t "$THREADS" \
      "$r1" "$r2" \
      > "$LOGDIR/kallisto_${sample}.stdout.log" 2> "$LOGDIR/kallisto_${sample}.stderr.log"
  fi
done

# ---------------------------
# 11) Summarize kallisto results to matrix (TPM & counts)
# ---------------------------
log "Step 8: Summarize kallisto outputs to matrices..."
TPM_MAT="$KALDIR/NR_TPM_matrix.tsv"
CNT_MAT="$KALDIR/NR_counts_matrix.tsv"

KALDIR="$KALDIR" python3 - <<'PY'
import os, glob
import pandas as pd

KALDIR=os.environ.get("KALDIR")
assert KALDIR and os.path.isdir(KALDIR), KALDIR

abunds=sorted(glob.glob(os.path.join(KALDIR, "*", "abundance.tsv")))
if not abunds:
    raise SystemExit("No abundance.tsv found under kallisto output directories.")

tpm = None
cnt = None

for f in abunds:
    sample=os.path.basename(os.path.dirname(f))
    df=pd.read_csv(f, sep="\t")
    df=df[["target_id","tpm","est_counts"]].copy().set_index("target_id")
    if tpm is None:
        tpm=df[["tpm"]].rename(columns={"tpm":sample})
        cnt=df[["est_counts"]].rename(columns={"est_counts":sample})
    else:
        tpm=tpm.join(df[["tpm"]].rename(columns={"tpm":sample}), how="outer")
        cnt=cnt.join(df[["est_counts"]].rename(columns={"est_counts":sample}), how="outer")

tpm=tpm.fillna(0.0)
cnt=cnt.fillna(0.0)

tpm_out=os.path.join(KALDIR, "NR_TPM_matrix.tsv")
cnt_out=os.path.join(KALDIR, "NR_counts_matrix.tsv")
tpm.to_csv(tpm_out, sep="\t")
cnt.to_csv(cnt_out, sep="\t")
print("[OK] wrote:", tpm_out)
print("[OK] wrote:", cnt_out)
PY

log "All done."
log "Key outputs:"
log "  NR genes:      $NR_GENES_FNA"
log "  NR proteins:   $NR_PROTS_FAA"
log "  eggNOG ann:    $EGGDIR/NR.emapper.annotations"
log "  kallisto TPM:  $TPM_MAT"
log "  kallisto cnt:  $CNT_MAT"
