#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

########################################
# Config (reuse existing genotype assets)
########################################
OUTDIR="/home/gao/data/host_microbe/hostgenome/GWAS_AP"

# New phenotype (MAG TPM table)
PHENO_TSV="/home/gao/data/host_microbe/hostgenome/phenotype/GWAS_MAG_TPM.txt"

# Optional extra covariates (if you have any besides PCs)
# Must contain: FID IID ...
COVAR_TSV="${COVAR_TSV:-}"   # leave empty by default

# CPU controls
THREADS="${THREADS:-128}"

# Genotype assets (already exist from your AP pipeline)
QC_PREFIX="$OUTDIR/03_qc/host_autosome.qc"
PCA_OUT="$OUTDIR/04_pca/host_autosome.pca"

# GWAS behavior
FORCE="${FORCE:-0}"                 # 0/1
PHENO_INVNRM="${PHENO_INVNRM:-1}"   # 1: use --pheno-quantile-normalize ; 0: raw TPM
N_PCS="${N_PCS:-10}"                # use PC1..PC10

########################################
# Helpers
########################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }
need plink2
need python3
need awk
need sort

have_plink_prefix() {
  local pref="$1"
  [[ -s "${pref}.pgen" && -s "${pref}.pvar" && -s "${pref}.psam" ]]
}

########################################
# Checks
########################################
echo "[check] OUTDIR    = $OUTDIR"
echo "[check] PHENO_TSV = $PHENO_TSV"
echo "[check] QC_PREFIX = $QC_PREFIX"
echo "[check] PCA_OUT   = $PCA_OUT"

[[ -s "$PHENO_TSV" ]] || { echo "ERROR: phenotype not found: $PHENO_TSV" >&2; exit 1; }
have_plink_prefix "$QC_PREFIX" || { echo "ERROR: missing QC pfile: ${QC_PREFIX}.pgen/.pvar/.psam" >&2; exit 1; }
[[ -s "${PCA_OUT}.eigenvec" ]] || { echo "ERROR: missing PCA eigenvec: ${PCA_OUT}.eigenvec" >&2; exit 1; }

mkdir -p "$OUTDIR"/{05_gwas_mag,06_clump_mag,07_annot_mag,00_lists_mag,log_mag}

########################################
# 1) Build KEEP list = intersect(pheno, genotype psam)
########################################
echo "[1] Build keep list (pheno ∩ genotype) ..."

KEEP="$OUTDIR/00_lists_mag/keep.MAG.txt"
PSAM="${QC_PREFIX}.psam"

if [[ -s "$KEEP" && "$FORCE" -eq 0 ]]; then
  echo "[1] Skip keep (exists): $KEEP"
else
  PHENO="$PHENO_TSV" PSAM="$PSAM" OUT="$KEEP" python3 - <<'PY'
import pandas as pd, os, sys
pheno = os.environ["PHENO"]
psam  = os.environ["PSAM"]
out   = os.environ["OUT"]

p = pd.read_csv(pheno, sep=r"\s+|\t", engine="python")
if not set(["FID","IID"]).issubset(p.columns):
    raise SystemExit("PHENO_TSV must contain columns: FID IID ...")

s = pd.read_csv(psam, sep=r"\s+|\t", engine="python")
# plink2 psam usually has #FID IID
s = s.rename(columns={"#FID":"FID"})
fid_col = [c for c in s.columns if c.lstrip("#").upper()=="FID"][0]
iid_col = [c for c in s.columns if c.upper()=="IID"][0]
s = s.rename(columns={fid_col:"FID", iid_col:"IID"})

m = p[["FID","IID"]].merge(s[["FID","IID"]], on=["FID","IID"], how="inner").drop_duplicates()
m.to_csv(out, sep="\t", index=False, header=False)
print("KEEP n =", m.shape[0], "->", out, file=sys.stderr)

if m.shape[0] < 50:
    print("WARNING: very small intersection. Check sample IDs.", file=sys.stderr)
PY
fi

########################################
# 2) Build covariate table: PC1..PCn (+ optional extra covars)
########################################
echo "[2] Build covariate table (PCs + optional covars) ..."

MERGED_COV="$OUTDIR/05_gwas_mag/covars.MAG.merged.tsv"

if [[ -s "$MERGED_COV" && "$FORCE" -eq 0 ]]; then
  echo "[2] Skip merged covars (exists): $MERGED_COV"
else
  KEEP_PATH="$KEEP" \
  PCA_PATH="${PCA_OUT}.eigenvec" \
  COV_PATH="${COVAR_TSV:-}" \
  OUT_PATH="$MERGED_COV" \
  N_PCS="$N_PCS" \
  python3 - <<'PY'
import pandas as pd, os, sys
keep_path = os.environ["KEEP_PATH"]
pca_path  = os.environ["PCA_PATH"]
cov_path  = os.environ.get("COV_PATH","")
out_path  = os.environ["OUT_PATH"]
n_pcs     = int(os.environ.get("N_PCS","10"))

keep = pd.read_csv(keep_path, sep="\t", header=None, names=["FID","IID"])

pca = pd.read_csv(pca_path, sep=r"\s+|\t", engine="python", header=None)
cols = ["FID","IID"] + [f"PC{i}" for i in range(1, pca.shape[1]-1)]
pca.columns = cols
pca = pca.merge(keep, on=["FID","IID"], how="inner")

pc_cols = ["FID","IID"] + [c for c in pca.columns if c.startswith("PC")][:n_pcs]
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
print("Merged covar n =", m.shape[0], "->", out_path, file=sys.stderr)
PY
fi

########################################
# 3) Extract phenotype column names (all MAG columns)
########################################
echo "[3] Detect MAG phenotype columns ..."

PHENO_NAMES_TXT="$OUTDIR/00_lists_mag/pheno.mag.names.txt"
PHENO_NAMES_CSV="$OUTDIR/00_lists_mag/pheno.mag.names.csv"

if [[ -s "$PHENO_NAMES_TXT" && -s "$PHENO_NAMES_CSV" && "$FORCE" -eq 0 ]]; then
  echo "[3] Skip pheno names (exists)."
else
  PHENO="$PHENO_TSV" OUT1="$PHENO_NAMES_TXT" OUT2="$PHENO_NAMES_CSV" python3 - <<'PY'
import pandas as pd, os
pheno = os.environ["PHENO"]
out1  = os.environ["OUT1"]
out2  = os.environ["OUT2"]

df = pd.read_csv(pheno, sep=r"\s+|\t", engine="python", nrows=5)
cols = [c for c in df.columns if c not in ("FID","IID")]
# 你这个文件是 MAG0637, MAG0657... 这种列名
with open(out1,"w") as w:
    for c in cols:
        w.write(c+"\n")
with open(out2,"w") as w:
    w.write(",".join(cols))
print("MAG phenotypes =", len(cols))
PY
fi

N_MAG=$(wc -l < "$PHENO_NAMES_TXT" | awk '{print $1}')
echo "[3] MAG traits detected: $N_MAG"
if [[ "$N_MAG" -lt 1 ]]; then
  echo "ERROR: No phenotype columns found besides FID/IID in $PHENO_TSV" >&2
  exit 1
fi

PHENO_NAME_LIST=$(cat "$PHENO_NAMES_CSV")

########################################
# 4) GWAS: Host SNPs ~ MAG TPM (continuous)
########################################
echo "[4] Run GWAS for all MAG traits ..."

GWAS_OUT="$OUTDIR/05_gwas_mag/MAG_gwas"

# build covar-name PC list
COVAR_NAME_LIST=$(python3 - <<PY
n=int("$N_PCS")
print(",".join([f"PC{i}" for i in range(1,n+1)]))
PY
)

GWAS_LOG="${GWAS_OUT}.log"
if [[ -s "$GWAS_LOG" && "$FORCE" -eq 0 ]]; then
  echo "[4] GWAS log exists; you can FORCE=1 to rerun: $GWAS_LOG"
else
  # phenotype normalization flag
  PHENO_NORM_FLAG=()
  if [[ "$PHENO_INVNRM" -eq 1 ]]; then
    PHENO_NORM_FLAG=(--pheno-quantile-normalize)
  fi

  plink2 \
    --pfile "$QC_PREFIX" \
    --chr-set 29 no-xy \
    --keep "$KEEP" \
    --pheno "$PHENO_TSV" \
    --pheno-name "$PHENO_NAME_LIST" \
    "${PHENO_NORM_FLAG[@]}" \
    --covar "$MERGED_COV" \
    --covar-name "$COVAR_NAME_LIST" \
    --covar-variance-standardize \
    --glm hide-covar \
    --out "$GWAS_OUT" \
    --threads "$THREADS"
fi

echo "[4] Done. Outputs in: $OUTDIR/05_gwas_mag/"
echo "    Example: ${GWAS_OUT}.MAG0637.glm.linear"

########################################
# 5) Optional: summarize significant hits per MAG (P<=1e-5)
########################################
echo "[5] Summarize suggestive hits per MAG ..."

SIG_SUM="$OUTDIR/05_gwas_mag/MAG_gwas.sig_P1e5.summary.tsv"
if [[ -s "$SIG_SUM" && "$FORCE" -eq 0 ]]; then
  echo "[5] Skip summary (exists): $SIG_SUM"
else
  GWAS_PREFIX="$GWAS_OUT" PHENOLIST="$PHENO_NAMES_TXT" OUT="$SIG_SUM" python3 - <<'PY'
import pandas as pd, os, glob

pref = os.environ["GWAS_PREFIX"]
phenolist = os.environ["PHENOLIST"]
out = os.environ["OUT"]

traits = [x.strip() for x in open(phenolist) if x.strip()]
rows = []
for t in traits:
    f = f"{pref}.{t}.glm.linear"
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    if "TEST" in df.columns:
        df = df[df["TEST"]=="ADD"]
    if "P" not in df.columns:
        continue
    n1 = int((df["P"]<=1e-5).sum())
    n2 = int((df["P"]<=1e-6).sum())
    minp = df["P"].min()
    rows.append((t, n1, n2, minp))

res = pd.DataFrame(rows, columns=["Trait","n_P<=1e-5","n_P<=1e-6","minP"]).sort_values("minP")
res.to_csv(out, sep="\t", index=False)
print("Wrote", out, "nTraits=", res.shape[0])
PY
fi

echo "All done!"
