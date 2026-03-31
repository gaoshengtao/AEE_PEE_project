#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

########################################
# Config
########################################
OUTDIR="/home/gao/data/host_microbe/hostgenome/GWAS_AP"

# existing genotype / PCA / phenotype
QC_PREFIX="$OUTDIR/03_qc/host_autosome.qc"
PCA_EIGENVEC="$OUTDIR/04_pca/host_autosome.pca.eigenvec"
PHENO_TSV="/home/gao/data/host_microbe/hostgenome/phenotype/GWAS_MAG_TPM.txt"

# existing MAG GWAS summary
GWAS_SIG_SUM="$OUTDIR/05_gwas_mag/MAG_gwas.sig_P1e5.summary.tsv"

# optional extra covariates besides PCs
COVAR_TSV="${COVAR_TSV:-}"   # can be empty

# outputs
H2_DIR="$OUTDIR/08_h2_mag"
LIST_DIR="$OUTDIR/00_lists_mag"

# params
N_PCS="${N_PCS:-10}"
THREADS="${THREADS:-32}"
FORCE="${FORCE:-0}"

# gcta executable
GCTA="/home/gao/apps/gcta-1.95.1-linux-x86_64/gcta64"

########################################
# Helpers
########################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }

need plink2
need python3
need awk
need sort
[[ -x "$GCTA" ]] || { echo "ERROR: GCTA not executable: $GCTA" >&2; exit 1; }

mkdir -p "$H2_DIR" "$LIST_DIR"

########################################
# Checks
########################################
[[ -s "${QC_PREFIX}.pgen" && -s "${QC_PREFIX}.pvar" && -s "${QC_PREFIX}.psam" ]] || {
  echo "ERROR: missing QC pfile: ${QC_PREFIX}.pgen/.pvar/.psam" >&2
  exit 1
}
[[ -s "$PCA_EIGENVEC" ]] || {
  echo "ERROR: missing PCA eigenvec: $PCA_EIGENVEC" >&2
  exit 1
}
[[ -s "$PHENO_TSV" ]] || {
  echo "ERROR: missing phenotype: $PHENO_TSV" >&2
  exit 1
}

########################################
# 1) Build GRM once
########################################
echo "[1] Build GRM ..."

BED_PREFIX="$H2_DIR/host_autosome.qc.for_gcta"
GRM_PREFIX="$H2_DIR/host_autosome.qc"

# 1.1 pfile -> bed/bim/fam (autosomes 1..29 only)
if [[ ! -s "${BED_PREFIX}.bed" || ! -s "${BED_PREFIX}.bim" || ! -s "${BED_PREFIX}.fam" || "$FORCE" -eq 1 ]]; then
  plink2 \
    --pfile "$QC_PREFIX" \
    --chr-set 29 no-xy \
    --chr 1-29 \
    --make-bed \
    --out "$BED_PREFIX" \
    --threads "$THREADS"
else
  echo "[1] Skip make-bed (exists): ${BED_PREFIX}.bed/.bim/.fam"
fi

# 1.2 make GRM
if [[ ! -s "${GRM_PREFIX}.grm.bin" || ! -s "${GRM_PREFIX}.grm.N.bin" || ! -s "${GRM_PREFIX}.grm.id" || "$FORCE" -eq 1 ]]; then
  "$GCTA" \
    --bfile "$BED_PREFIX" \
    --autosome-num 29 \
    --make-grm \
    --thread-num "$THREADS" \
    --out "$GRM_PREFIX"
else
  echo "[1] Skip make-grm (exists): ${GRM_PREFIX}.grm.*"
fi

########################################
# 2) Detect MAG phenotype columns
########################################
echo "[2] Detect MAG phenotype columns ..."

PHENO_NAMES_TXT="$LIST_DIR/pheno.mag.names.txt"

if [[ ! -s "$PHENO_NAMES_TXT" || "$FORCE" -eq 1 ]]; then
  PHENO="$PHENO_TSV" OUT1="$PHENO_NAMES_TXT" python3 - <<'PY'
import pandas as pd, os
pheno = os.environ["PHENO"]
out1 = os.environ["OUT1"]

df = pd.read_csv(pheno, sep=r"\s+|\t", engine="python", nrows=5)
cols = [c for c in df.columns if c not in ("FID", "IID")]
with open(out1, "w") as w:
    for c in cols:
        w.write(c + "\n")
print("MAG phenotypes =", len(cols))
PY
else
  echo "[2] Skip detect phenotype columns (exists): $PHENO_NAMES_TXT"
fi

########################################
# 3) Build keep list and qcovar
########################################
echo "[3] Build qcovar file (PCs + optional numeric covars) ..."

Q_COVAR="$H2_DIR/MAG.qcovar.tsv"
KEEP="$LIST_DIR/keep.MAG.txt"
PSAM="${QC_PREFIX}.psam"

# 3.1 keep list = intersection(pheno, psam)
if [[ ! -s "$KEEP" || "$FORCE" -eq 1 ]]; then
  PHENO="$PHENO_TSV" PSAM="$PSAM" OUT="$KEEP" python3 - <<'PY'
import pandas as pd, os
pheno = os.environ["PHENO"]
psam = os.environ["PSAM"]
out = os.environ["OUT"]

p = pd.read_csv(pheno, sep=r"\s+|\t", engine="python")
if not set(["FID","IID"]).issubset(p.columns):
    raise SystemExit("PHENO_TSV must contain columns: FID IID ...")

s = pd.read_csv(psam, sep=r"\s+|\t", engine="python")
s = s.rename(columns={"#FID":"FID"})
fid_col = [c for c in s.columns if c.lstrip("#").upper()=="FID"][0]
iid_col = [c for c in s.columns if c.upper()=="IID"][0]
s = s.rename(columns={fid_col:"FID", iid_col:"IID"})

m = p[["FID","IID"]].merge(s[["FID","IID"]], on=["FID","IID"], how="inner").drop_duplicates()
m.to_csv(out, sep="\t", index=False, header=False)
print("KEEP n =", m.shape[0])
PY
else
  echo "[3] Skip keep list (exists): $KEEP"
fi

# 3.2 qcovar
if [[ ! -s "$Q_COVAR" || "$FORCE" -eq 1 ]]; then
  KEEP_PATH="$KEEP" \
  PCA_PATH="$PCA_EIGENVEC" \
  COV_PATH="${COVAR_TSV:-}" \
  OUT_PATH="$Q_COVAR" \
  N_PCS="$N_PCS" \
  python3 - <<'PY'
import pandas as pd, os

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
pc_cols = ["FID","IID"] + [f"PC{i}" for i in range(1, n_pcs+1)]
pca = pca[pc_cols]

if cov_path and os.path.exists(cov_path) and os.path.getsize(cov_path) > 0:
    cov = pd.read_csv(cov_path, sep=r"\s+|\t", engine="python")
    if not set(["FID","IID"]).issubset(cov.columns):
        raise SystemExit("COVAR_TSV must contain columns: FID IID ...")
    cov = cov.merge(keep, on=["FID","IID"], how="inner")
    extra_cols = [c for c in cov.columns if c not in ("FID","IID")]
    num_cols = [c for c in extra_cols if pd.api.types.is_numeric_dtype(cov[c])]
    m = pca.merge(cov[["FID","IID"] + num_cols], on=["FID","IID"], how="left")
else:
    m = pca

m.to_csv(out_path, sep="\t", index=False, header=False)
print("Wrote qcovar:", out_path, "n=", m.shape[0], "cols=", m.shape[1])
PY
else
  echo "[3] Skip qcovar (exists): $Q_COVAR"
fi

########################################
# 4) Run REML for each MAG
########################################
echo "[4] Run SNP-heritability (REML) for each MAG ..."

while read -r trait; do
  [[ -n "$trait" ]] || continue

  PHENO_ONE="$H2_DIR/${trait}.pheno.tsv"
  OUT_ONE="$H2_DIR/${trait}"

  # 4.1 extract one phenotype
  if [[ ! -s "$PHENO_ONE" || "$FORCE" -eq 1 ]]; then
    TRAIT="$trait" PHENO="$PHENO_TSV" KEEP="$KEEP" OUT="$PHENO_ONE" python3 - <<'PY'
import pandas as pd, os

trait = os.environ["TRAIT"]
pheno = os.environ["PHENO"]
keepf = os.environ["KEEP"]
out = os.environ["OUT"]

df = pd.read_csv(pheno, sep=r"\s+|\t", engine="python")
keep = pd.read_csv(keepf, sep="\t", header=None, names=["FID","IID"])

if trait not in df.columns:
    raise SystemExit(f"Trait not found: {trait}")

x = df[["FID","IID",trait]].merge(keep, on=["FID","IID"], how="inner")
x = x.dropna(subset=[trait])
x.to_csv(out, sep="\t", index=False, header=False)
print("Wrote", out, "n=", x.shape[0])
PY
  fi

  # 4.2 REML
  if [[ -s "${OUT_ONE}.hsq" && "$FORCE" -eq 0 ]]; then
    echo "[4] Skip $trait (exists): ${OUT_ONE}.hsq"
  else
    echo "[4] REML: $trait"
    "$GCTA" \
      --grm "$GRM_PREFIX" \
      --pheno "$PHENO_ONE" \
      --qcovar "$Q_COVAR" \
      --autosome-num 29 \
      --reml \
      --reml-no-constrain \
      --thread-num "$THREADS" \
      --out "$OUT_ONE"
  fi

done < "$PHENO_NAMES_TXT"

########################################
# 5) Summarize h2 results
########################################
echo "[5] Summarize h2 results ..."

H2_SUM="$H2_DIR/MAG_h2_summary.tsv"

python3 - <<'PY'
import os, glob, pandas as pd, re

h2_dir = "/home/gao/data/host_microbe/hostgenome/GWAS_AP/08_h2_mag"
rows = []

for f in sorted(glob.glob(os.path.join(h2_dir, "*.hsq"))):
    trait = os.path.basename(f).replace(".hsq", "")

    vg = ve = vp = h2 = h2_se = lrt = pval = None

    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 2:
                continue

            key = parts[0]
            if key == "V(G)":
                vg = parts[1]
            elif key == "V(e)":
                ve = parts[1]
            elif key == "Vp":
                vp = parts[1]
            elif key == "V(G)/Vp":
                h2 = parts[1]
                if len(parts) > 2:
                    h2_se = parts[2]
            elif key == "LRT":
                lrt = parts[1]
            elif key == "Pval":
                pval = parts[1]

    rows.append((trait, vg, ve, vp, h2, h2_se, lrt, pval))

res = pd.DataFrame(rows, columns=[
    "Trait","V_G","V_e","V_p","h2_SNP","SE","LRT","Pval"
])

for c in ["V_G","V_e","V_p","h2_SNP","SE","LRT","Pval"]:
    res[c] = pd.to_numeric(res[c], errors="coerce")

res = res.sort_values(["h2_SNP","Pval"], ascending=[False, True])
out = os.path.join(h2_dir, "MAG_h2_summary.tsv")
res.to_csv(out, sep="\t", index=False)
print("Wrote", out, "n=", res.shape[0])
PY

########################################
# 6) Merge h2 with GWAS signal summary
########################################
echo "[6] Merge h2 with GWAS signal summary ..."

H2_GWAS_SUM="$H2_DIR/MAG_h2_vs_GWAS_signal.tsv"

if [[ -s "$GWAS_SIG_SUM" ]]; then
  H2_SUM_IN="$H2_SUM" GWAS_SUM_IN="$GWAS_SIG_SUM" OUT_SUM="$H2_GWAS_SUM" python3 - <<'PY'
import os
import numpy as np
import pandas as pd

h2_file = os.environ["H2_SUM_IN"]
gwas_file = os.environ["GWAS_SUM_IN"]
out_file = os.environ["OUT_SUM"]

h2 = pd.read_csv(h2_file, sep="\t")
gwas = pd.read_csv(gwas_file, sep="\t")

# 兼容列名
if "Trait" not in gwas.columns:
    raise SystemExit("GWAS summary missing column: Trait")
if "minP" not in gwas.columns:
    raise SystemExit("GWAS summary missing column: minP")

m = h2.merge(gwas, on="Trait", how="outer")

m["minus_log10_minP"] = np.where(
    (m["minP"].notna()) & (m["minP"] > 0),
    -np.log10(m["minP"]),
    np.nan
)

# 更适合后续分析的排序
m = m.sort_values(
    ["h2_SNP", "minus_log10_minP", "Trait"],
    ascending=[False, False, True]
)

m.to_csv(out_file, sep="\t", index=False)
print("Wrote", out_file, "n=", m.shape[0])
PY
else
  echo "[6] Skip merge: GWAS summary not found -> $GWAS_SIG_SUM"
fi

echo
echo "Done!"
echo "Heritability results:"
echo "  $H2_DIR"
echo "Summary:"
echo "  $H2_SUM"
if [[ -s "$H2_GWAS_SUM" ]]; then
  echo "h2 vs GWAS signal:"
  echo "  $H2_GWAS_SUM"
fi
