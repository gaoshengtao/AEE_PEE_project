#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(stringi)
})

counts_TSV  <- "/home/gao/data/host_microbe/gene_catalog/05_kallisto/NR_counts_matrix.tsv"
ANNO_TSV <- "/home/gao/data/host_microbe/gene_catalog/04_eggnog/NR.emapper.annotations"
OUTDIR   <- "/home/gao/data/host_microbe/gene_catalog/06_aggregate_counts_by_eggnog"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# 需要聚合的字段（每个字段输出一张矩阵）
FIELDS <- c(
  "eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC",
  "KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC",
  "CAZy","BiGG_Reaction","PFAMs"
)

# -------------------------
# helper: 规范化分隔符并拆分
# -------------------------
split_tokens <- function(x) {
  # 统一空白
  x <- gsub("\\s+", "", x)
  # eggNOG/KEGG 经常用逗号分隔；有时也可能有 ; | 作为分隔
  # 注意：eggNOG_OGs 里也有 '|'（层级信息），我们不按 '|' 拆，只按 , ; 分
  # PFAMs 通常用逗号分隔
  stri_split_regex(x, pattern = "[,;]+", omit_empty = TRUE, simplify = FALSE)
}

# -------------------------
# 1) 读注释（只读我们要用的列，省内存）
# -------------------------
anno_cols_need <- c("#query", FIELDS)
anno <- fread(
  ANNO_TSV,
  sep = "\t",
  header = TRUE,
  quote = "",
  data.table = TRUE,
  showProgress = TRUE,
  select = intersect(anno_cols_need, names(fread(ANNO_TSV, nrows = 0)))
)

if (!"#query" %in% names(anno)) stop("Annotation file does not have #query column.")
setnames(anno, "#query", "gene_id")

# 只保留 gene_id + 目标字段
anno <- anno[, c("gene_id", intersect(FIELDS, names(anno))), with = FALSE]

# -------------------------
# 2) 读 counts 矩阵（wide：gene × sample）
# -------------------------
counts_dt <- fread(counts_TSV, sep = "\t", header = TRUE, quote = "", data.table = TRUE, showProgress = TRUE)

# 找到第一列 gene_id（通常是 target_id / gene / #query 之类）
gene_col <- names(counts_dt)[1]
setnames(counts_dt, gene_col, "gene_id")

# 取样本列
sample_cols <- setdiff(names(counts_dt), "gene_id")
if (length(sample_cols) < 2) stop("counts matrix seems to have <2 sample columns.")

# 与注释取交集 gene（避免无注释的基因拖累）
common_gene <- intersect(counts_dt$gene_id, anno$gene_id)
cat("[INFO] genes in counts:", nrow(counts_dt), "\n")
cat("[INFO] genes in anno:", nrow(anno), "\n")
cat("[INFO] common genes :", length(common_gene), "\n")

counts_dt <- counts_dt[gene_id %in% common_gene]
anno   <- anno[gene_id %in% common_gene]

# 建立 gene 索引（用于稀疏矩阵乘法）
gene_levels <- counts_dt$gene_id
gene_index  <- seq_along(gene_levels)
names(gene_index) <- gene_levels

# 把 counts 变成 dense 矩阵：gene × sample（注意：很大，需足够内存）
counts_mat <- as.matrix(counts_dt[, ..sample_cols])
storage.mode(counts_mat) <- "double"

# -------------------------
# 3) 对每个字段：展开 + 分摊 + 稀疏聚合
# -------------------------
for (field in FIELDS) {
  if (!field %in% names(anno)) {
    cat("[WARN] field not found in anno:", field, " (skip)\n")
    next
  }
  
  cat("\n[INFO] Aggregating field:", field, "\n")
  
  dtf <- anno[, .(gene_id, val = get(field))]
  dtf <- dtf[!is.na(val) & val != "" & val != "-"]
  
  if (nrow(dtf) == 0) {
    cat("[WARN] no valid annotations for field:", field, "\n")
    next
  }
  
  # 拆分 token（每个 gene 多个 val）
  toks <- split_tokens(dtf$val)
  dt_long <- data.table(
    gene_id = rep(dtf$gene_id, lengths(toks)),
    term    = unlist(toks, use.names = FALSE)
  )
  dt_long <- dt_long[term != "" & term != "-"]
  
  # 每个 gene 在该字段下的 term 个数 k -> 权重 1/k
  k_dt <- dt_long[, .N, by = gene_id]
  setnames(k_dt, "N", "k")
  dt_long <- dt_long[k_dt, on = "gene_id"]
  dt_long[, w := 1 / k]
  
  # 去重（同一 gene-term 如果重复出现，保留一次）
  dt_long <- unique(dt_long, by = c("gene_id", "term"))
  
  # term 索引
  term_levels <- sort(unique(dt_long$term))
  term_id <- match(dt_long$term, term_levels)
  
  # gene 索引
  g_id <- gene_index[dt_long$gene_id]
  
  # 构建稀疏映射：term × gene 的权重矩阵
  # M[term, gene] = weight
  M <- sparseMatrix(
    i = term_id,
    j = g_id,
    x = dt_long$w,
    dims = c(length(term_levels), length(gene_levels))
  )
  
  # 聚合：term × sample = (term × gene) %*% (gene × sample)
  agg_mat <- as.matrix(M %*% counts_mat)
  colnames(agg_mat) <- sample_cols
  rownames(agg_mat) <- term_levels
  
  # 输出
  out_tsv <- file.path(OUTDIR, paste0("counts_by_", field, ".tsv"))
  out_dt <- data.table(term = rownames(agg_mat))
  out_dt <- cbind(out_dt, as.data.table(agg_mat))
  fwrite(out_dt, out_tsv, sep = "\t")
  
  # 输出一个简单统计
  cat("[INFO] field:", field, " terms:", nrow(agg_mat), " output:", out_tsv, "\n")
}

cat("\n[DONE] All aggregations written to:", OUTDIR, "\n")

