#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

########################################
# Config
########################################
BAM_BASE="/home/gao/data/host_microbe/hostgenome/bwa_out"
REF_FA="/home/gao/backup/database/kneaddata/cow/Bos_taurus.ARS-UCD1.3.dna.toplevel.fasta"
OUTDIR="/home/gao/data/host_microbe/hostgenome/GWAS_AP"

# GWAS inputs
PHENO_TSV="/home/gao/data/host_microbe/hostgenome/phenotype/AP_pheno.tsv"  # FID IID AP
COVAR_TSV="/home/gao/data/host_microbe/hostgenome/phenotype/AP_covar.tsv"  # optional FID IID ...

# CPU controls
THREADS="${THREADS:-128}"
MAX_JOBS="${MAX_JOBS:-29}"
MP_THREADS="${MP_THREADS:-2}"

# mpileup controls
MAX_DEPTH="${MAX_DEPTH:-800}"
MIN_MQ="${MIN_MQ:-30}"
MIN_BQ="${MIN_BQ:-20}"

# Genotype/site filters
MIN_DP="${MIN_DP:-3}"
MIN_QUAL="${MIN_QUAL:-30}"

# PLINK QC
MAF="${MAF:-0.01}"
GENO="${GENO:-0.1}"
HWE="${HWE:-1e-6}"
MIND="${MIND:-0.1}"

# Force rerun per-chr calling
FORCE="${FORCE:-0}"   # 0/1

########################################
# Helpers
########################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }
need samtools
need bcftools
need plink2
need python3
need awk
need find
need sort
need bedtools
need gzip


mkdir -p "$OUTDIR"/{00_lists,01_vcf,02_plink,03_qc,04_pca,05_gwas,log}

# skip helper: check existence of a "prefix" trio
have_plink_prefix() {
  local pref="$1"
  [[ -s "${pref}.pgen" && -s "${pref}.pvar" && -s "${pref}.psam" ]]
}

########################################
# Reference checks
########################################
if [[ ! -s "$REF_FA" ]]; then
  echo "ERROR: REF_FA not found: $REF_FA" >&2
  exit 1
fi
if [[ ! -s "${REF_FA}.fai" ]]; then
  echo "[ref] samtools faidx $REF_FA"
  samtools faidx "$REF_FA"
fi

########################################
# 0) Collect ALL BAMs
########################################
echo "[0] Collect ALL BAMs under: $BAM_BASE"
BAM_LIST="$OUTDIR/00_lists/bams.all.list"
find "$BAM_BASE" -mindepth 2 -maxdepth 2 -type f -name "*.sorted.bam" | sort > "$BAM_LIST"

N_BAMS=$(wc -l < "$BAM_LIST" | awk '{print $1}')
echo "BAMs found: $N_BAMS"
if [[ "$N_BAMS" -lt 100 ]]; then
  echo "ERROR: Too few BAMs found. Check BAM_BASE structure." >&2
  exit 1
fi

ALL_IDS="$OUTDIR/00_lists/sample_all_bam.txt"
awk -F'/' '{print $(NF-1)}' "$BAM_LIST" | sort -u > "$ALL_IDS"
echo "Unique sample IDs from BAMs: $(wc -l < "$ALL_IDS")"

########################################
# 1) Ensure BAM index
########################################
echo "[1] Ensure BAM index ..."
while IFS= read -r bam; do
  [[ -s "${bam}.bai" || -s "${bam%.bam}.bai" ]] || samtools index -@ "$THREADS" "$bam"
done < "$BAM_LIST"

########################################
# 2) Build autosome chr list (1..29 only)
########################################
echo "[2-prep] Build autosome chr list (1..29 only) ..."
CHR_LIST_ALL="$OUTDIR/00_lists/chr.all.fai.list"
cut -f1 "${REF_FA}.fai" > "$CHR_LIST_ALL"

CHR_LIST="$OUTDIR/00_lists/chr.autosome.1_29.list"
awk '
function inrange(x){ return (x>=1 && x<=29) }
{
  chr=$1
  if (chr ~ /^[0-9]+$/) {
    n=chr+0
    if (inrange(n)) print chr
  } else if (chr ~ /^chr[0-9]+$/) {
    sub(/^chr/,"",chr)
    n=chr+0
    if (inrange(n)) print $1
  }
}
' "$CHR_LIST_ALL" > "$CHR_LIST"

N_CHR=$(wc -l < "$CHR_LIST" | awk '{print $1}')
echo "Autosome contigs selected: $N_CHR"
if [[ "$N_CHR" -lt 10 ]]; then
  echo "ERROR: autosome list too short. REF naming unexpected; check ${REF_FA}.fai" >&2
  head -n 20 "${REF_FA}.fai" >&2
  exit 1
fi

########################################
# 2) Joint variant calling by chromosome
########################################
echo "[2] Joint calling by chromosome (autosomes only) ..."
CHR_DIR="$OUTDIR/01_vcf/by_chr"
mkdir -p "$CHR_DIR"

RAW_VCF="$OUTDIR/01_vcf/raw.autosome.bcf"

run_one_chr() {
  set -Eeuo pipefail
  local chr="$1"
  local out_bcf="$CHR_DIR/chr.${chr}.bcf"
  local out_csi="${out_bcf}.csi"

  if [[ "$FORCE" -eq 0 && -s "$out_bcf" && -s "$out_csi" ]]; then
    echo "[chr] skip $chr (exists)"
    return 0
  fi

  echo "[chr] start $chr"
  bcftools mpileup \
    -r "$chr" \
    -f "$REF_FA" \
    -b "$BAM_LIST" \
    -Ou \
    -a FORMAT/DP,FORMAT/AD \
    --max-depth "$MAX_DEPTH" \
    --threads "$MP_THREADS" \
    -q "$MIN_MQ" \
    -Q "$MIN_BQ" \
  | bcftools call -m -Ob -o "$out_bcf"

  bcftools index -f "$out_bcf"
  echo "[chr] done  $chr"
}
export -f run_one_chr
export REF_FA BAM_LIST CHR_DIR MAX_DEPTH MP_THREADS MIN_MQ MIN_BQ FORCE

jobcnt=0
while read -r chr; do
  run_one_chr "$chr" &
  jobcnt=$((jobcnt+1))
  if [[ "$jobcnt" -ge "$MAX_JOBS" ]]; then
    wait -n
    jobcnt=$((jobcnt-1))
  fi
done < "$CHR_LIST"
wait

########################################
# 2b) Concatenate per-chr BCFs -> RAW_VCF (SKIP if exists)
########################################
if [[ -s "$RAW_VCF" && -s "${RAW_VCF}.csi" && "$FORCE" -eq 0 ]]; then
  echo "[2] Skip concat (exists): $RAW_VCF"
else
  echo "[2] Concatenate per-chr BCFs -> $RAW_VCF"
  BCF_LIST="$OUTDIR/00_lists/bychr_bcf.autosome.list"
  : > "$BCF_LIST"
  while read -r chr; do
    echo "$CHR_DIR/chr.${chr}.bcf" >> "$BCF_LIST"
  done < "$CHR_LIST"

  missing=0
  while read -r f; do
    [[ -s "$f" ]] || { echo "ERROR: missing per-chr bcf: $f" >&2; missing=1; }
  done < "$BCF_LIST"
  [[ "$missing" -eq 0 ]] || exit 1

  bcftools concat -f "$BCF_LIST" -Ob -o "$RAW_VCF"
  bcftools index -f "$RAW_VCF"
fi

########################################
# 3) Filter for GWAS (SKIP if exists)
########################################
FLT_VCF="$OUTDIR/01_vcf/filtered.autosome.bcf"
if [[ -s "$FLT_VCF" && -s "${FLT_VCF}.csi" && "$FORCE" -eq 0 ]]; then
  echo "[3] Skip filtering (exists): $FLT_VCF"
else
  echo "[3] GWAS-friendly filtering ..."
  bcftools view -m2 -M2 -v snps "$RAW_VCF" -Ou \
  | bcftools filter -Ou -S . -e "FMT/DP<${MIN_DP}" \
  | bcftools filter -Ou -e "QUAL<${MIN_QUAL}" \
  | bcftools view -Ob -o "$FLT_VCF"

  bcftools index -f "$FLT_VCF"
fi

echo "[3-check] Header DP fields:"
bcftools view -h "$FLT_VCF" | grep -E '##(FORMAT|INFO)=<ID=DP' -n || true

############################################
# 4.0 生成样本映射表
################################################
SAMPLE_RENAME="$OUTDIR/00_lists/sample.rename.tsv"

python3 - <<'PY'
import pandas as pd

bam_list = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/00_lists/bams.all.list"
out = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/00_lists/sample.rename.tsv"

rows=[]
with open(bam_list) as f:
    for line in f:
        bam=line.strip()
        if not bam:
            continue
        sid=bam.split("/")[-2]  # G0619
        # oldFID oldIID newFID newIID
        rows.append((bam, bam, sid, sid))

pd.DataFrame(rows).drop_duplicates().to_csv(out, sep="\t", header=False, index=False)
print("Wrote", out, "n=", len(rows))
PY

########################################
# 4) VCF -> PLINK2 (SKIP if exists)
#    FIX: enforce unique variant IDs so --indep-pairwise won't fail
########################################
echo "[4] VCF -> PLINK2 ..."
PLINK_PREFIX="$OUTDIR/02_plink/host_autosome"

if have_plink_prefix "$PLINK_PREFIX" && [[ "$FORCE" -eq 0 ]]; then
  echo "[4] Skip PLINK convert (exists): ${PLINK_PREFIX}.pgen/.pvar/.psam"
else
  plink2 \
    --bcf "$FLT_VCF" \
    --chr-set 29 no-xy \
    --double-id \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 500 \
    --update-ids "$SAMPLE_RENAME" \
    --rm-dup force-first \
    --make-pgen \
    --out "$PLINK_PREFIX" \
    --threads "$THREADS"
fi

########################################
# 5) QC in PLINK2 (SKIP if exists)
########################################
echo "[5] QC ..."
QC_PREFIX="$OUTDIR/03_qc/host_autosome.qc"

if have_plink_prefix "$QC_PREFIX" && [[ "$FORCE" -eq 0 ]]; then
  echo "[5] Skip QC (exists): ${QC_PREFIX}.pgen/.pvar/.psam"
else
  plink2 \
    --pfile "$PLINK_PREFIX" \
    --chr-set 29 no-xy \
    --mind "$MIND" \
    --geno "$GENO" \
    --maf "$MAF" \
    --hwe "$HWE" \
    --make-pgen \
    --out "$QC_PREFIX" \
    --threads "$THREADS"
fi

########################################
# 6) LD pruning + PCA
########################################
echo "[6] PCA ..."
PRUNE_OUT="$OUTDIR/04_pca/host_autosome.prune"
PCA_OUT="$OUTDIR/04_pca/host_autosome.pca"

# prune (skip if exists)
if [[ -s "${PRUNE_OUT}.prune.in" && -s "${PRUNE_OUT}.prune.out" && "$FORCE" -eq 0 ]]; then
  echo "[6] Skip pruning (exists): ${PRUNE_OUT}.prune.in"
else
  plink2 \
    --pfile "$QC_PREFIX" \
    --chr-set 29 no-xy \
    --indep-pairwise 200 50 0.2 \
    --out "$PRUNE_OUT" \
    --threads "$THREADS"
fi

# pca (skip if exists)
if [[ -s "${PCA_OUT}.eigenvec" && -s "${PCA_OUT}.eigenval" && "$FORCE" -eq 0 ]]; then
  echo "[6] Skip PCA (exists): ${PCA_OUT}.eigenvec"
else
  plink2 \
    --pfile "$QC_PREFIX" \
    --chr-set 29 no-xy \
    --extract "${PRUNE_OUT}.prune.in" \
    --pca approx 20 \
    --out "$PCA_OUT" \
    --threads "$THREADS"
fi

echo "PCA written: ${PCA_OUT}.eigenvec"

########################################
# 7) GWAS: AP phenotype (INTERSECT at GWAS stage)
########################################
echo "[7] GWAS (AP) using intersect of phenotype & genotype ..."

KEEP="$OUTDIR/00_lists/keep.AP.txt"
PHENO="$PHENO_TSV" \
PSAM="${QC_PREFIX}.psam" \
OUT="$KEEP" \
python3 - <<'PY'
import pandas as pd, os
pheno = os.environ["PHENO"]
psam  = os.environ["PSAM"]
out   = os.environ["OUT"]
p = pd.read_csv(pheno, sep=r"\s+|\t", engine="python")
s = pd.read_csv(psam,  sep=r"\s+|\t", engine="python")
s = s.rename(columns={"#FID":"FID"})
fid_col = [c for c in s.columns if c.lstrip("#").upper()=="FID"][0]
iid_col = [c for c in s.columns if c.upper()=="IID"][0]
s = s.rename(columns={fid_col:"FID", iid_col:"IID"})
m = p[["FID","IID"]].merge(s[["FID","IID"]], on=["FID","IID"], how="inner").drop_duplicates()
m.to_csv(out, sep="\t", index=False, header=False)
print("KEEP n =", m.shape[0], "->", out)
PY

MERGED_COV="$OUTDIR/05_gwas/covars.AP.merged.tsv"
KEEP_PATH="$KEEP" \
PCA_PATH="${PCA_OUT}.eigenvec" \
COV_PATH="$COVAR_TSV" \
OUT_PATH="$MERGED_COV" \
python3 - <<'PY'
import pandas as pd, os
keep_path = os.environ["KEEP_PATH"]
pca_path  = os.environ["PCA_PATH"]
cov_path  = os.environ.get("COV_PATH","")
out_path  = os.environ["OUT_PATH"]

keep = pd.read_csv(keep_path, sep="\t", header=None, names=["FID","IID"])

pca = pd.read_csv(pca_path, sep=r"\s+|\t", engine="python", header=None)
cols = ["FID","IID"] + [f"PC{i}" for i in range(1, pca.shape[1]-1)]
pca.columns = cols
pca = pca.merge(keep, on=["FID","IID"], how="inner")
pc_cols = ["FID","IID"] + [c for c in pca.columns if c.startswith("PC")][:10]
pca = pca[pc_cols]

if cov_path and os.path.exists(cov_path) and os.path.getsize(cov_path) > 0:
    cov = pd.read_csv(cov_path, sep=r"\s+|\t", engine="python")
    if not set(["FID","IID"]).issubset(cov.columns):
        raise SystemExit("COVAR_TSV must contain columns: FID IID ...")
    cov = cov.merge(keep, on=["FID","IID"], how="inner")
    m = pca.merge(cov, on=["FID","IID"], how="left")
else:
    m = pca

m.to_csv(out_path, sep="\t", index=False)
print("Merged covar n =", m.shape[0], "->", out_path)
PY

GWAS_OUT="$OUTDIR/05_gwas/AP_gwas"
plink2 \
  --pfile "$QC_PREFIX" \
  --chr-set 29 no-xy \
  --keep "$KEEP" \
  --pheno "$PHENO_TSV" \
  --pheno-name AP \
  --covar "$MERGED_COV" \
  --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --glm hide-covar \
  --covar-variance-standardize \
  --out "$GWAS_OUT" \
  --threads "$THREADS"
  
  
  
########################################
# 8) LD clumping (SKIP if exists)
########################################
echo "[8] LD clumping ..."

GWAS_GLM="${GWAS_OUT}.AP.glm.linear"
CLUMP_DIR="$OUTDIR/06_clump"
mkdir -p "$CLUMP_DIR"

# 你的 GWAS 输出列第一列是 #CHROM，先统一去掉 '#'
GWAS_FOR_CLUMP="$CLUMP_DIR/AP_gwas.for_clump.tsv"
CLUMP_OUT="$CLUMP_DIR/AP_gwas.clump"

# 生成 clump 输入文件（plink2 --clump 需要至少 SNP + P）
if [[ -s "$GWAS_FOR_CLUMP" && "$FORCE" -eq 0 ]]; then
  echo "[8] Skip prepare clump input (exists): $GWAS_FOR_CLUMP"
else
  echo "[8] Prepare clump input: $GWAS_FOR_CLUMP"
  python3 - <<'PY'
import pandas as pd

infile = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/05_gwas/AP_gwas.AP.glm.linear"
outfile = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/06_clump/AP_gwas.for_clump.tsv"

df = pd.read_csv(infile, sep="\t")
# rename #CHROM -> CHROM if present
df.columns = [c.lstrip("#") for c in df.columns]

# 只要 ADD
if "TEST" in df.columns:
    df = df[df["TEST"]=="ADD"]

# 必须有 ID 和 P
df = df[["ID","P"]].dropna()
df.to_csv(outfile, sep="\t", index=False)
print("Wrote", outfile, "n=", df.shape[0])
PY
fi

# clumping
if [[ -s "${CLUMP_OUT}.clumped" && "$FORCE" -eq 0 ]]; then
  echo "[8] Skip clumping (exists): ${CLUMP_OUT}.clumped"
else
  plink2 \
    --pfile "$QC_PREFIX" \
    --chr-set 29 no-xy \
    --clump "$GWAS_FOR_CLUMP" \
    --clump-snp-field ID \
    --clump-field P \
    --clump-p1 1e-5 \
    --clump-p2 1e-3 \
    --clump-r2 0.2 \
    --clump-kb 500 \
    --out "$CLUMP_OUT" \
    --threads "$THREADS"
fi

echo "[8] Clump results: ${CLUMP_OUT}.clumped"


########################################
# 9) Annotate Top SNPs (P<=1e-5) with nearest gene using GTF (SKIP if exists)
########################################
echo "[9] Annotate Top SNPs (P<=1e-5) -> nearest gene ..."

GTF="/home/gao/backup/database/kneaddata/cow/Bos_taurus.ARS-UCD1.3.113.chr.gtf"
ANN_DIR="$OUTDIR/07_annot"
mkdir -p "$ANN_DIR"

TOP_SNP_TSV="$ANN_DIR/AP.top1e5.snps.tsv"
TOP_SNP_BED="$ANN_DIR/AP.top1e5.snps.bed"
GENE_BED="$ANN_DIR/bos_taurus.genes.bed"
ANN_OUT="$ANN_DIR/AP.top1e5.nearest_gene.tsv"

# 9.1 提取 top SNP（P<=1e-5）
if [[ -s "$TOP_SNP_TSV" && "$FORCE" -eq 0 ]]; then
  echo "[9] Skip top SNP extract (exists): $TOP_SNP_TSV"
else
  python3 - <<'PY'
import pandas as pd

infile = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/05_gwas/AP_gwas.AP.glm.linear"
outfile = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot/AP.top1e5.snps.tsv"

df = pd.read_csv(infile, sep="\t")
df.columns = [c.lstrip("#") for c in df.columns]

# ADD only
if "TEST" in df.columns:
    df = df[df["TEST"]=="ADD"]

# filter
df = df[df["P"].notna()]
df = df[df["P"] <= 1e-5].copy()

# keep useful cols if exist
keep_cols = [c for c in ["CHROM","POS","ID","REF","ALT","A1","A1_FREQ","BETA","SE","T_STAT","P","OBS_CT"] if c in df.columns]
df = df[keep_cols].sort_values("P")

df.to_csv(outfile, sep="\t", index=False)
print("Top1e-5 n =", df.shape[0], "->", outfile)
PY
fi

# 9.2 SNP TSV -> BED (0-based)
if [[ -s "$TOP_SNP_BED" && "$FORCE" -eq 0 ]]; then
  echo "[9] Skip SNP BED (exists): $TOP_SNP_BED"
else
  awk 'BEGIN{OFS="\t"} NR>1{chrom=$1; pos=$2; id=$3; start=pos-1; end=pos; print chrom,start,end,id}' \
    "$TOP_SNP_TSV" > "$TOP_SNP_BED"
fi

# ✅ 关键：无论是否 skip，都把 SNP BED 排序一次（bedtools 要求）
sort -k1,1 -k2,2n "$TOP_SNP_BED" -o "$TOP_SNP_BED"


# 9.3 从 GTF 提取 gene 区间 -> BED
# 只取 feature == "gene"，提 gene_id 和 gene_name（如果有）
if [[ -s "$GENE_BED" && "$FORCE" -eq 0 ]]; then
  echo "[9] Skip gene BED (exists): $GENE_BED"
else
  python3 - <<'PY'
import re, sys

gtf = "/home/gao/backup/database/kneaddata/cow/Bos_taurus.ARS-UCD1.3.113.chr.gtf"
out = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot/bos_taurus.genes.bed"

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
        # bed: chrom start end name strand gene_id
        w.write("\t".join([chrom, str(start0), str(end1), gene_name, strand, gene_id]) + "\n")

print("Wrote", out, file=sys.stderr)
PY
fi

# ✅ 关键：无论是否 skip，都把 gene BED 排序一次
sort -k1,1 -k2,2n "$GENE_BED" -o "$GENE_BED"


# 9.4 最近基因注释（带距离）
# 输出字段：snp_chr snp_start snp_end snp_id | gene_chr gene_start gene_end gene_name strand gene_id | distance
if [[ -s "$ANN_OUT" && "$FORCE" -eq 0 ]]; then
  echo "[9] Skip annotation (exists): $ANN_OUT"
else
  bedtools closest -a "$TOP_SNP_BED" -b "$GENE_BED" -d > "$ANN_DIR/tmp.closest.tsv"
  

  # 加表头，并与 GWAS 信息合并（把 P/BETA 等拼进去）
  python3 - <<'PY'
import pandas as pd

top = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot/AP.top1e5.snps.tsv"
closest = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot/tmp.closest.tsv"
out = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/07_annot/AP.top1e5.nearest_gene.tsv"

g = pd.read_csv(top, sep="\t")
# columns already CHROM POS ID ...
g = g.rename(columns={"CHROM":"snp_chr","POS":"snp_pos","ID":"snp_id"})

c = pd.read_csv(closest, sep="\t", header=None)
c.columns = [
    "snp_chr","snp_start","snp_end","snp_id",
    "gene_chr","gene_start","gene_end","gene_name","gene_strand","gene_id",
    "dist"
]

# 合并
m = g.merge(c[["snp_id","gene_name","gene_id","gene_strand","gene_start","gene_end","dist"]],
            on="snp_id", how="left")

m.to_csv(out, sep="\t", index=False)
print("Annotated ->", out, "n=", m.shape[0])
PY

  rm -f "$ANN_DIR/tmp.closest.tsv"
fi

echo "[9] Annotation output: $ANN_OUT"
echo "[9] Done."


echo "Done!"
echo "Genotype assets (autosomes only):"
echo "  RAW BCF : $RAW_VCF"
echo "  FLT BCF : $FLT_VCF"
echo "  PLINK   : ${QC_PREFIX}.pgen/.pvar/.psam"
echo "  PCA     : ${PCA_OUT}.eigenvec"
echo "AP GWAS result:"
echo "  ${GWAS_OUT}.AP.glm.linear"

