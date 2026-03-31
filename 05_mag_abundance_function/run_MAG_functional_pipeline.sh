#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

############################################
# Input / output
############################################
REP_DIR="/home/gao/data/host_microbe/MAGs/drep/dereplicated_genomes"
GENOME_INFO="/home/gao/data/host_microbe/MAGs/drep/data_tables/genomeInformation.csv"

OUTDIR="/home/gao/data/host_microbe/MAGs/function_annotation_eggnog"
THREADS="${THREADS:-64}"
JOBS="${JOBS:-6}"          # 并行 MAG 数（JOBS*EMAPPER_CPU 不要超过CPU）
PROD_MODE="${PROD_MODE:-single}"  # single / meta

############################################
# eggNOG-mapper (your paths)
############################################
EGGNOG_MAPPER="$HOME/apps/eggnog-mapper-2.1.13/emapper.py"
EGGNOG_DB_DIR="/home/gao/databases/eggnog"
EGGNOG_CONDA_ENV="${EGGNOG_CONDA_ENV:-eggnog}"
EMAPPER_CPU="${EMAPPER_CPU:-8}"

############################################
# Tools
############################################
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }
need prodigal
need python3
need awk
need find
need sort
need xargs
need conda

# emapper path existence check (no PATH needed)
[[ -s "$EGGNOG_MAPPER" ]] || { echo "ERROR: emapper.py not found: $EGGNOG_MAPPER" >&2; exit 1; }
[[ -d "$EGGNOG_DB_DIR" ]] || { echo "ERROR: eggNOG DB dir not found: $EGGNOG_DB_DIR" >&2; exit 1; }

mkdir -p "$OUTDIR"/{00_lists,01_prodigal,02_eggnog,03_matrices,log,tmp}

############################################
# 0) List repMAG FASTA files
############################################
FASTALIST="$OUTDIR/00_lists/repMAGs.list"
find "$REP_DIR" -maxdepth 1 -type f \( -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" \) \
  | sort > "$FASTALIST"

N=$(wc -l < "$FASTALIST" | awk '{print $1}')
echo "[0] repMAG FASTA files: $N"
[[ "$N" -ge 10 ]] || { echo "ERROR: too few MAGs in $REP_DIR" >&2; exit 1; }

############################################
# 1) Prodigal per MAG
############################################
echo "[1] Prodigal gene calling ..."

prodigal_one(){
  set -Eeuo pipefail
  local fa="$1"
  local mag
  mag="$(basename "$fa")"; mag="${mag%.*}"

  local out="$OUTDIR/01_prodigal/$mag"
  mkdir -p "$out"

  local aa="$out/${mag}.faa"
  local nucl="$out/${mag}.fna"
  local gff="$out/${mag}.gff"

  if [[ -s "$aa" && -s "$gff" ]]; then
    echo "[SKIP] prodigal $mag"
    return 0
  fi

  prodigal -i "$fa" -a "$aa" -d "$nucl" -f gff -o "$gff" -p "$PROD_MODE" \
    2> "$OUTDIR/log/${mag}.prodigal.log"

  echo "[DONE] prodigal $mag"
}
export -f prodigal_one

# ✅ 关键：把 OUTDIR/PROD_MODE 导出给子进程
export OUTDIR PROD_MODE

cat "$FASTALIST" | xargs -n 1 -P "$JOBS" bash -c 'prodigal_one "$@"' _

############################################
# 2) eggNOG-mapper per MAG (conda run)
############################################
echo "[2] eggNOG-mapper annotation ..."

eggnog_one(){
  set -Eeuo pipefail
  local fa="$1"
  local mag
  mag="$(basename "$fa")"; mag="${mag%.*}"

  local aa="$OUTDIR/01_prodigal/$mag/${mag}.faa"
  local outdir="$OUTDIR/02_eggnog/$mag"
  mkdir -p "$outdir"

  # eggNOG-mapper 2.1.x 默认输出：<out>.emapper.annotations(.tsv)
  local ann1="$outdir/${mag}.emapper.annotations"
  local ann2="$outdir/${mag}.emapper.annotations.tsv"

  if [[ -s "$ann1" || -s "$ann2" ]]; then
    echo "[SKIP] eggnog $mag"
    return 0
  fi
  [[ -s "$aa" ]] || { echo "[WARN] missing proteins: $aa" >&2; return 0; }

  echo "[RUN] eggnog $mag"
    conda run -n "$EGGNOG_CONDA_ENV" \
    python "$EGGNOG_MAPPER" \
      -i "$aa" \
      --itype proteins \
      --output "$mag" \
      --output_dir "$outdir" \
      --data_dir "$EGGNOG_DB_DIR" \
      --cpu "$EMAPPER_CPU" \
      --override \
      > "$OUTDIR/log/${mag}.eggnog.log" 2>&1


  echo "[DONE] eggnog $mag"
}
export -f eggnog_one

# ✅ 关键：把这些变量导出给子进程
export OUTDIR EGGNOG_MAPPER EGGNOG_DB_DIR EGGNOG_CONDA_ENV EMAPPER_CPU

cat "$FASTALIST" | xargs -n 1 -P "$JOBS" bash -c 'eggnog_one "$@"' _

############################################
# 3) Build MAG x function matrices (KO / EC / COGcat / eggNOG_OGs)
############################################
echo "[3] Build MAG x function matrices ..."

# ✅ 关键：把 GENOME_INFO 真正传给 python
OUTDIR="$OUTDIR" GENOME_INFO="$GENOME_INFO" python3 - <<'PY'
import os, glob, re, io
import pandas as pd

OUTDIR = os.environ["OUTDIR"]
EGGDIR = os.path.join(OUTDIR, "02_eggnog")
META_IN = os.environ.get("GENOME_INFO","")

os.makedirs(os.path.join(OUTDIR, "03_matrices"), exist_ok=True)

# ---------------------------
# Read eggNOG annotations:
#  - drop lines starting with "##"
#  - keep header line "#query\t..."
#  - remove leading "#" in header
# ---------------------------
def read_emapper_annotations(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = []
        for ln in f:
            if ln.startswith("##"):
                continue
            lines.append(ln)
    if not lines:
        return pd.DataFrame()

    if lines[0].startswith("#"):
        lines[0] = lines[0].lstrip("#")  # "#query" -> "query"

    buf = io.StringIO("".join(lines))
    return pd.read_csv(buf, sep="\t", low_memory=False)

# -------- genome size from repMAG fasta --------
fastalist = os.path.join(OUTDIR, "00_lists", "repMAGs.list")
sizes = {}
with open(fastalist) as f:
    for line in f:
        fa = line.strip()
        if not fa:
            continue
        mag = os.path.splitext(os.path.basename(fa))[0]
        total = 0
        with open(fa, "r", encoding="utf-8", errors="replace") as fh:
            for s in fh:
                if s.startswith(">"):
                    continue
                total += len(s.strip())
        sizes[mag] = total

meta = pd.DataFrame({"MAG": list(sizes.keys()),
                     "genome_size_bp": list(sizes.values())})
meta["genome_size_Mb"] = meta["genome_size_bp"]/1e6

# -------- merge genomeInformation.csv if provided --------
if META_IN and os.path.exists(META_IN) and os.path.getsize(META_IN)>0:
    gi = pd.read_csv(META_IN)
    cand = [c for c in gi.columns if c.lower() in ["genome","genome_id","bin_id","name","mag","id"]]
    idcol = cand[0] if cand else gi.columns[0]
    gi = gi.rename(columns={idcol:"MAG"})
    gi["MAG"] = gi["MAG"].astype(str).str.replace(r"\.(fa|fna|fasta|fas|fsa)$", "", regex=True, case=False)
    meta = meta.merge(gi, on="MAG", how="left")

meta.to_csv(os.path.join(OUTDIR, "03_matrices", "MAG_metadata.tsv"),
            sep="\t", index=False)

ko_pat = re.compile(r"K\d{5}$")

ko_records, ec_records, cog_records, og_records = [], [], [], []
pathway_records, module_records = [], []   # ✅ 新增

def split_multi(raw: str):
    # 兼容逗号/分号/竖线分隔（eggNOG 有时会用逗号，有时会是逗号+空格）
    raw = str(raw)
    for sep in [";", "|"]:
        raw = raw.replace(sep, ",")
    return [x.strip() for x in raw.split(",") if x.strip()]

def normalize_kegg_list_items(items):
    """
    eggNOG-mapper 的 KEGG_Pathway/Module 有时会长这样：
      path:map00010, map00020
      ko00010
      M00001
    这里做一个轻量清洗：去掉前缀 'path:' 'map' 'ko' 等，仅保留核心ID
    - Pathway: 保留 mapXXXXX 或 koXXXXX（以防返回 koXXXXX）
    - Module:  保留 MXXXXX
    """
    out = []
    for it in items:
        it = it.strip()
        if not it or it in ["-", "--"]:
            continue
        it = it.replace("path:", "").replace("PATH:", "")
        out.append(it)
    return out

# scan eggnog outputs
for magdir in sorted(glob.glob(os.path.join(EGGDIR, "*"))):
    if not os.path.isdir(magdir):
        continue
    mag = os.path.basename(magdir)

    ann = None
    for fn in [f"{mag}.emapper.annotations", f"{mag}.emapper.annotations.tsv"]:
        p = os.path.join(magdir, fn)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            ann = p
            break
    if not ann:
        continue

    df = read_emapper_annotations(ann)
    if df.empty:
        continue

    # columns are: query, ... , eggNOG_OGs, COG_category, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, ...
    colmap = {c.lower(): c for c in df.columns}
    ko_col  = colmap.get("kegg_ko")
    ec_col  = colmap.get("ec")
    cog_col = colmap.get("cog_category")
    og_col  = colmap.get("eggnog_ogs")

    # ✅ 新增：Pathway / Module（大小写兼容）
    path_col = colmap.get("kegg_pathway") or colmap.get("kegg_pathways")
    modu_col = colmap.get("kegg_module")  or colmap.get("kegg_modules")

    if ko_col is None and ec_col is None and cog_col is None and og_col is None and path_col is None and modu_col is None:
        print(f"[WARN] No KO/EC/COG/OG/Pathway/Module columns in {mag}. First cols:", list(df.columns)[:15])
        continue

    for _, r in df.iterrows():
        # --- KO ---
        if ko_col and pd.notna(r.get(ko_col)):
            raw = str(r[ko_col]).strip()
            if raw and raw not in ["-", "--"]:
                raw = raw.replace("ko:", "")
                for ko in split_multi(raw):
                    if ko_pat.match(ko):
                        ko_records.append((mag, ko))

        # --- EC ---
        if ec_col and pd.notna(r.get(ec_col)):
            raw = str(r[ec_col]).strip()
            if raw and raw not in ["-", "--"]:
                for ec in split_multi(raw):
                    if ec and ec not in ["-", "--"]:
                        ec_records.append((mag, ec))

        # --- COG category letters ---
        if cog_col and pd.notna(r.get(cog_col)):
            s = str(r[cog_col]).strip()
            if s and s not in ["-", "--"]:
                for ch in s:
                    if ch.isalpha():
                        cog_records.append((mag, ch))

        # --- eggNOG OGs (fine-grained) ---
        if og_col and pd.notna(r.get(og_col)):
            raw = str(r[og_col]).strip()
            if raw and raw not in ["-", "--"]:
                # 例子：COG0001@1,ENOG4105C20@2759
                for item in split_multi(raw):
                    if item in ["-", "--"]:
                        continue
                    og = item.split("@", 1)[0].strip()
                    if og and og not in ["-", "--"]:
                        og_records.append((mag, og))

        # ✅ --- KEGG Pathway ---
        if path_col and pd.notna(r.get(path_col)):
            raw = str(r[path_col]).strip()
            if raw and raw not in ["-", "--"]:
                items = normalize_kegg_list_items(split_multi(raw))
                for pw in items:
                    if pw and pw not in ["-", "--"]:
                        pathway_records.append((mag, pw))

        # ✅ --- KEGG Module ---
        if modu_col and pd.notna(r.get(modu_col)):
            raw = str(r[modu_col]).strip()
            if raw and raw not in ["-", "--"]:
                items = normalize_kegg_list_items(split_multi(raw))
                for md in items:
                    if md and md not in ["-", "--"]:
                        module_records.append((mag, md))

def write_mat(records, colname, outname):
    if not records:
        print("EMPTY:", outname)
        return
    d = pd.DataFrame(records, columns=["MAG", colname])
    mat = d.value_counts(["MAG", colname]).reset_index(name="copy_number")
    wide = mat.pivot(index="MAG", columns=colname, values="copy_number").fillna(0).astype(int)
    wide.to_csv(outname, sep="\t")
    print("WROTE:", outname, wide.shape)

write_mat(ko_records,       "KO",         os.path.join(OUTDIR, "03_matrices", "MAGxKO_copy_number.tsv"))
write_mat(ec_records,       "EC",         os.path.join(OUTDIR, "03_matrices", "MAGxEC_copy_number.tsv"))
write_mat(cog_records,      "COG",        os.path.join(OUTDIR, "03_matrices", "MAGxCOGcat_copy_number.tsv"))
write_mat(og_records,       "eggNOG_OG",   os.path.join(OUTDIR, "03_matrices", "MAGxeggNOG_OGs_copy_number.tsv"))
# ✅ 新增输出
write_mat(pathway_records,  "KEGG_Pathway",os.path.join(OUTDIR, "03_matrices", "MAGxKEGG_Pathway_copy_number.tsv"))
write_mat(module_records,   "KEGG_Module", os.path.join(OUTDIR, "03_matrices", "MAGxKEGG_Module_copy_number.tsv"))
PY


echo "[DONE] Outputs:"
echo "  $OUTDIR/03_matrices/MAG_metadata.tsv"
echo "  $OUTDIR/03_matrices/MAGxKO_copy_number.tsv"
echo "  $OUTDIR/03_matrices/MAGxEC_copy_number.tsv"
echo "  $OUTDIR/03_matrices/MAGxCOGcat_copy_number.tsv"
echo "  $OUTDIR/03_matrices/MAGxeggNOG_OGs_copy_number.tsv"


