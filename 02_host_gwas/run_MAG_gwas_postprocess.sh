#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

########################################
# Config
########################################
OUTDIR="/home/gao/data/host_microbe/hostgenome/GWAS_AP"

# genotype assets from finished AP pipeline
QC_PREFIX="$OUTDIR/03_qc/host_autosome.qc"

# existing MAG GWAS results
GWAS_DIR="$OUTDIR/05_gwas_mag"
GWAS_PREFIX="$GWAS_DIR/MAG_gwas"

# annotation
GTF="/home/gao/backup/database/kneaddata/cow/Bos_taurus.ARS-UCD1.3.113.chr.gtf"

# outputs
CLUMP_DIR="$OUTDIR/06_clump_mag"
ANN_DIR="$OUTDIR/07_annot_mag"
LIST_DIR="$OUTDIR/00_lists_mag"

# clump params
CLUMP_P1="${CLUMP_P1:-1e-5}"
CLUMP_P2="${CLUMP_P2:-1e-3}"
CLUMP_R2="${CLUMP_R2:-0.2}"
CLUMP_KB="${CLUMP_KB:-500}"

# threads / behavior
THREADS="${THREADS:-32}"
FORCE="${FORCE:-0}"

########################################
# Helpers
########################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }
need plink2
need python3
need bedtools
need sort
need awk

mkdir -p "$CLUMP_DIR" "$ANN_DIR" "$LIST_DIR"

########################################
# Checks
########################################
[[ -s "${QC_PREFIX}.pgen" && -s "${QC_PREFIX}.pvar" && -s "${QC_PREFIX}.psam" ]] || {
  echo "ERROR: missing QC pfile: ${QC_PREFIX}.pgen/.pvar/.psam" >&2
  exit 1
}

[[ -s "$GTF" ]] || {
  echo "ERROR: missing GTF: $GTF" >&2
  exit 1
}

GWAS_FILES=( "$GWAS_DIR"/MAG_gwas.*.glm.linear )
if [[ ${#GWAS_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no MAG GWAS files found under $GWAS_DIR" >&2
  exit 1
fi

echo "[check] MAG GWAS files found: ${#GWAS_FILES[@]}"

########################################
# 6) LD clumping for each MAG
########################################
echo "[6] LD clumping for each MAG ..."

for f in "${GWAS_FILES[@]}"; do
  base=$(basename "$f")
  trait="${base#MAG_gwas.}"
  trait="${trait%.glm.linear}"

  CLUMP_IN="$CLUMP_DIR/${trait}.for_clump.tsv"
  CLUMP_OUT="$CLUMP_DIR/${trait}"

  # prepare clump input
  if [[ -s "$CLUMP_IN" && "$FORCE" -eq 0 ]]; then
    echo "[6] Skip clump input: $CLUMP_IN"
  else
    echo "[6] Prepare clump input for $trait"
    INFILE="$f" OUTFILE="$CLUMP_IN" python3 - <<'PY'
import os
import pandas as pd

infile = os.environ["INFILE"]
outfile = os.environ["OUTFILE"]

df = pd.read_csv(infile, sep="\t")
df.columns = [c.lstrip("#") for c in df.columns]

if "TEST" in df.columns:
    df = df[df["TEST"] == "ADD"]

need = ["ID", "P"]
miss = [x for x in need if x not in df.columns]
if miss:
    raise SystemExit(f"Missing columns in {infile}: {', '.join(miss)}")

df = df[["ID", "P"]].dropna()
df = df[(df["P"] > 0) & (df["P"] <= 1)]
df.to_csv(outfile, sep="\t", index=False)
print(f"Wrote {outfile} n={df.shape[0]}")
PY
  fi

  # clump
  if [[ -s "${CLUMP_OUT}.clumps" && "$FORCE" -eq 0 ]]; then
    echo "[6] Skip clumping for $trait (exists): ${CLUMP_OUT}.clumps"
  else
    echo "[6] Clumping $trait ..."
    plink2 \
      --pfile "$QC_PREFIX" \
      --chr-set 29 no-xy \
      --clump "$CLUMP_IN" \
      --clump-snp-field ID \
      --clump-field P \
      --clump-p1 "$CLUMP_P1" \
      --clump-p2 "$CLUMP_P2" \
      --clump-r2 "$CLUMP_R2" \
      --clump-kb "$CLUMP_KB" \
      --out "$CLUMP_OUT" \
      --threads "$THREADS"
  fi
done

echo "[6] Done."

########################################
# 6b) Summarize clump results
########################################
echo "[6b] Summarize clump results ..."

CLUMP_SUM="$CLUMP_DIR/MAG_gwas.clump.summary.tsv"

if [[ -s "$CLUMP_SUM" && "$FORCE" -eq 0 ]]; then
  echo "[6b] Skip clump summary: $CLUMP_SUM"
else
  python3 - <<'PY'
import os, glob, pandas as pd

clump_dir = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/06_clump_mag"
rows = []

for f in sorted(glob.glob(os.path.join(clump_dir, "*.clumps"))):
    trait = os.path.basename(f).replace(".clumps", "")
    try:
        df = pd.read_csv(f, sep=r"\s+", engine="python", comment=None)
    except Exception:
        rows.append((trait, 0, None))
        continue

    if df.shape[0] == 0:
        rows.append((trait, 0, None))
        continue

    df.columns = [c.lstrip("#") for c in df.columns]
    lead_id = df.loc[0, "ID"] if "ID" in df.columns else None
    rows.append((trait, df.shape[0], lead_id))

res = pd.DataFrame(rows, columns=["Trait", "n_lead_loci", "top_lead_snp"])
res = res.sort_values(["n_lead_loci", "Trait"], ascending=[False, True])
out = os.path.join(clump_dir, "MAG_gwas.clump.summary.tsv")
res.to_csv(out, sep="\t", index=False)
print("Wrote", out, "n=", res.shape[0])
PY
fi

########################################
# 7) Prepare shared gene BED from GTF
########################################
echo "[7] Prepare gene BED from GTF ..."

GENE_BED="$ANN_DIR/bos_taurus.genes.bed"

if [[ -s "$GENE_BED" && "$FORCE" -eq 0 ]]; then
  echo "[7] Skip gene BED: $GENE_BED"
else
  GTF_IN="$GTF" BED_OUT="$GENE_BED" python3 - <<'PY'
import os, re

gtf = os.environ["GTF_IN"]
out = os.environ["BED_OUT"]

def get_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else ""

with open(gtf) as f, open(out, "w") as w:
    for line in f:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attr = parts
        if feature != "gene":
            continue
        start0 = int(start) - 1
        end1 = int(end)
        gene_id = get_attr(attr, "gene_id")
        gene_name = get_attr(attr, "gene_name")
        if not gene_name:
            gene_name = gene_id
        w.write("\t".join([chrom, str(start0), str(end1), gene_name, strand, gene_id]) + "\n")
print("Wrote", out)
PY
fi

sort -k1,1 -k2,2n "$GENE_BED" -o "$GENE_BED"

########################################
# 7b) Annotate top SNPs (P<=1e-5) for each MAG
########################################
echo "[7b] Annotate nearest genes for each MAG ..."

for f in "${GWAS_FILES[@]}"; do
  base=$(basename "$f")
  trait="${base#MAG_gwas.}"
  trait="${trait%.glm.linear}"

  TOP_SNP_TSV="$ANN_DIR/${trait}.top1e5.snps.tsv"
  TOP_SNP_BED="$ANN_DIR/${trait}.top1e5.snps.bed"
  ANN_OUT="$ANN_DIR/${trait}.top1e5.nearest_gene.tsv"

  # extract top SNPs
  if [[ -s "$TOP_SNP_TSV" && "$FORCE" -eq 0 ]]; then
    echo "[7b] Skip top SNP extract for $trait"
  else
    echo "[7b] Extract top SNPs for $trait"
    INFILE="$f" OUTFILE="$TOP_SNP_TSV" python3 - <<'PY'
import os
import pandas as pd

infile = os.environ["INFILE"]
outfile = os.environ["OUTFILE"]

df = pd.read_csv(infile, sep="\t")
df.columns = [c.lstrip("#") for c in df.columns]

if "TEST" in df.columns:
    df = df[df["TEST"] == "ADD"]

df = df[df["P"].notna()]
df = df[(df["P"] <= 1e-5) & (df["P"] > 0)].copy()

keep_cols = [c for c in [
    "CHROM","POS","ID","REF","ALT","A1","A1_FREQ","BETA","SE","T_STAT","P","OBS_CT"
] if c in df.columns]

df = df[keep_cols].sort_values("P")
df.to_csv(outfile, sep="\t", index=False)
print("Wrote", outfile, "n=", df.shape[0])
PY
  fi

  # if no top SNPs, skip annotation
  if [[ ! -s "$TOP_SNP_TSV" ]] || [[ $(wc -l < "$TOP_SNP_TSV") -le 1 ]]; then
    echo "[7b] No SNPs with P<=1e-5 for $trait, skip annotation."
    continue
  fi

  # tsv -> bed
  if [[ -s "$TOP_SNP_BED" && "$FORCE" -eq 0 ]]; then
    echo "[7b] Skip SNP BED for $trait"
  else
    awk 'BEGIN{OFS="\t"} NR>1{chrom=$1; pos=$2; id=$3; start=pos-1; end=pos; print chrom,start,end,id}' \
      "$TOP_SNP_TSV" > "$TOP_SNP_BED"
  fi

  sort -k1,1 -k2,2n "$TOP_SNP_BED" -o "$TOP_SNP_BED"

  # closest annotation
  if [[ -s "$ANN_OUT" && "$FORCE" -eq 0 ]]; then
    echo "[7b] Skip annotation for $trait"
  else
    echo "[7b] Annotate $trait"
    TMP_CLOSEST="$ANN_DIR/${trait}.tmp.closest.tsv"
    bedtools closest -a "$TOP_SNP_BED" -b "$GENE_BED" -d > "$TMP_CLOSEST"

    TOPFILE="$TOP_SNP_TSV" CLOSESTFILE="$TMP_CLOSEST" OUTFILE="$ANN_OUT" python3 - <<'PY'
import os
import pandas as pd

top = os.environ["TOPFILE"]
closest = os.environ["CLOSESTFILE"]
out = os.environ["OUTFILE"]

g = pd.read_csv(top, sep="\t")
g = g.rename(columns={"CHROM":"snp_chr","POS":"snp_pos","ID":"snp_id"})

c = pd.read_csv(closest, sep="\t", header=None)
c.columns = [
    "snp_chr","snp_start","snp_end","snp_id",
    "gene_chr","gene_start","gene_end","gene_name","gene_strand","gene_id",
    "dist"
]

m = g.merge(
    c[["snp_id","gene_name","gene_id","gene_strand","gene_start","gene_end","dist"]],
    on="snp_id",
    how="left"
)

m.to_csv(out, sep="\t", index=False)
print("Annotated ->", out, "n=", m.shape[0])
PY

    rm -f "$TMP_CLOSEST"
  fi
done

########################################
# 7c) Merge annotation summaries
########################################
echo "[7c] Merge all annotation summaries ..."

ANN_MERGED="$ANN_DIR/MAG_gwas.top1e5.nearest_gene.merged.tsv"

if [[ -s "$ANN_MERGED" && "$FORCE" -eq 0 ]]; then
  echo "[7c] Skip merged annotation: $ANN_MERGED"
else
  python3 - <<'PY'
import os, glob, pandas as pd

ann_dir = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot_mag"
files = sorted(glob.glob(os.path.join(ann_dir, "*.top1e5.nearest_gene.tsv")))
dfs = []
for f in files:
    trait = os.path.basename(f).replace(".top1e5.nearest_gene.tsv", "")
    try:
        df = pd.read_csv(f, sep="\t")
    except Exception:
        continue
    if df.shape[0] == 0:
        continue
    df.insert(0, "Trait", trait)
    dfs.append(df)

out = os.path.join(ann_dir, "MAG_gwas.top1e5.nearest_gene.merged.tsv")
if dfs:
    res = pd.concat(dfs, axis=0, ignore_index=True)
    res.to_csv(out, sep="\t", index=False)
    print("Wrote", out, "n=", res.shape[0])
else:
    print("No annotation tables found.")
PY
fi

echo
echo "Done!"
echo "Clump results:"
echo "  $CLUMP_DIR"
echo "Annotation results:"
echo "  $ANN_DIR"
