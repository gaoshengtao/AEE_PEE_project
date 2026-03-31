#!/usr/bin/env Rscript

## =========================================================
## MAG functional clustering by OG/KO copy-number (microbiome-style)
## - completeness correction only (no genome size correction)
## - compositional normalization: CLR (recommended) or CSS (optional)
## - feature filtering: prevalence + variance
## - dimension reduction: PCA + (optional) PCoA Bray-Curtis
## - clustering: HDBSCAN (no preset K) + DBSCAN optional
## - visualization: PCA/UMAP/tSNE/PCoA
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(uwot)
  library(Rtsne)
  library(vegan)
  library(cluster)
  library(dbscan)
  library(compositions)  # CLR
})

set.seed(42)

## -----------------------------
## User config
## -----------------------------
IN_FEATURE  <- "./微生物/MAGxKO_copy_number.tsv"  # 或 MAGxKO_copy_number.tsv
IN_META     <- "微生物/genomeInformation.csv"             # 至少包含 MAG + completeness
OUTDIR      <- "微生物/MAG_COG_cluster"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "Figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "Tables"),  showWarnings = FALSE, recursive = TRUE)

USE_METHOD  <- "CLR"  # "CLR" (推荐) 或 "CSS"（可选，需要 metagenomeSeq）
PSEUDOCOUNT <- 0.5    # CLR 需要 pseudocount，0.5 或 1 都行

# feature filtering
MIN_PREV_FRAC <- 0.05   # 至少在 5% 的 MAG 中出现（>0）
MIN_VAR_Q     <- 0.25   # 过滤低方差：保留 variance 位于上四分之三（可改 0.1/0.2/0.3）

# PCA dims
N_PCS <- 20

# UMAP / tSNE
UMAP_NN   <- 30
UMAP_MD   <- 0.2
TSNE_PERP <- 30

# HDBSCAN (no preset K)
# minPts 越大越“保守”，越容易把小簇当噪声；你可以试 10/15/20/30
HDBSCAN_MINPTS <- 20

# DBSCAN optional (if you want)
RUN_DBSCAN <- TRUE

## -----------------------------
## Helpers
## -----------------------------
msg <- function(...) cat(sprintf(...), "\n")

strip_ext <- function(x){
  sub("\\.(fa|fna|fasta|fas|fsa)$", "", x, ignore.case = TRUE)
}

read_meta_completeness <- function(meta_path){
  meta <- fread(meta_path, data.table = FALSE)
  # 容错列名
  mag_candidates  <- c("MAG","mag","genome","genome_id","Genome","ID","id")
  comp_candidates <- c("completeness","Completeness","checkm_completeness","completion")
  
  mag_col  <- mag_candidates[mag_candidates %in% colnames(meta)][1]
  comp_col <- comp_candidates[comp_candidates %in% colnames(meta)][1]
  if (is.na(mag_col))  stop("meta 里找不到 MAG 列（MAG/genome_id 等）")
  if (is.na(comp_col)) stop("meta 里找不到 completeness 列（completeness/Completeness 等）")
  
  meta <- meta %>%
    rename(MAG = all_of(mag_col)) %>%
    mutate(MAG = as.character(MAG),
           MAG = strip_ext(MAG),
           completeness_raw = as.numeric(.data[[comp_col]])) %>%
    mutate(comp_frac = ifelse(completeness_raw > 1.5, completeness_raw/100, completeness_raw)) %>%
    filter(!is.na(comp_frac), comp_frac > 0)
  
  meta %>% select(MAG, comp_frac) %>% distinct()
}

# CLR on rows (per MAG composition)
clr_rows <- function(mat, pseudocount = 0.5){
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0
  mat <- mat + pseudocount
  # compositions::clr expects compositions object
  # we want CLR across features per sample (row)
  # compositions::clr works on columns by default; transpose twice
  t(clr(acomp(t(mat))))
}

# CSS optional (metagenomeSeq)
css_rows <- function(mat){
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("USE_METHOD='CSS' 需要 metagenomeSeq：install.packages('BiocManager'); BiocManager::install('metagenomeSeq')")
  }
  suppressPackageStartupMessages(library(metagenomeSeq))
  # metagenomeSeq expects features x samples
  m <- t(as.matrix(mat))
  pheno <- AnnotatedDataFrame(data.frame(row.names = colnames(m)))
  feat  <- AnnotatedDataFrame(data.frame(row.names = rownames(m)))
  obj <- newMRexperiment(m, phenoData = pheno, featureData = feat)
  obj <- cumNorm(obj, p = cumNormStatFast(obj))
  norm <- MRcounts(obj, norm = TRUE, log = FALSE) # features x samples
  t(norm)
}

# silhouette helper (for any clustering; noise = 0 will be removed)
silhouette_score <- function(labels, dist_obj){
  lab <- as.integer(as.factor(labels))
  # silhouette requires >=2 clusters
  if (length(unique(lab)) < 2) return(NA_real_)
  mean(silhouette(lab, dist_obj)[, "sil_width"])
}

## -----------------------------
## 1) Read feature matrix
## -----------------------------
feat <- fread(IN_FEATURE, data.table = FALSE)
colnames(feat)[1] <- "MAG"
feat$MAG <- strip_ext(as.character(feat$MAG))

# numeric feature columns
feat_cols <- setdiff(colnames(feat), "MAG")
feat <- feat %>% mutate(across(all_of(feat_cols), ~ suppressWarnings(as.numeric(.x))))
feat[is.na(feat)] <- 0

meta <- read_meta_completeness(IN_META)
dat  <- feat %>% inner_join(meta, by = "MAG")
msg("[1] MAGs after merge = %d", nrow(dat))


# Build matrix
X_raw <- dat %>%
  select(MAG, all_of(feat_cols)) %>%
  tibble::column_to_rownames("MAG") %>%
  as.matrix()

comp <- dat$comp_frac
names(comp) <- dat$MAG

## -----------------------------
## 2) Completeness correction (only)
## copy_corrected = copy / comp_frac
## -----------------------------
X_corr <- X_raw
X_corr <- sweep(X_corr, 1, comp[rownames(X_corr)], `/`)
X_corr[!is.finite(X_corr)] <- 0

fwrite(data.frame(MAG=rownames(X_corr), X_corr, check.names = FALSE),
       file.path(OUTDIR, "Tables", "MAGxFeature_copy_corrected.tsv"),
       sep = "\t")

## -----------------------------
## 3) Feature filtering: prevalence + variance
## -----------------------------
prev <- colMeans(X_corr > 0)
keep_prev <- prev >= MIN_PREV_FRAC

vars <- apply(X_corr, 2, var)
cut_var <- quantile(vars[is.finite(vars)], probs = MIN_VAR_Q, na.rm = TRUE)
keep_var <- vars >= cut_var

keep <- keep_prev & keep_var
X_filt <- X_corr[, keep, drop = FALSE]

msg("[3] Features: raw=%d, after filter=%d (prev>=%.2f, var>=Q%.2f)",
    ncol(X_corr), ncol(X_filt), MIN_PREV_FRAC, MIN_VAR_Q)

fwrite(data.frame(feature=colnames(X_corr), prevalence=prev, variance=vars,
                  keep=keep),
       file.path(OUTDIR, "Tables", "feature_filtering_summary.tsv"),
       sep="\t")

## -----------------------------
## 4) Normalization (CLR recommended; CSS optional)
## -----------------------------
if (USE_METHOD == "CLR") {
  msg("[4] Normalization = CLR (pseudocount=%.3f)", PSEUDOCOUNT)
  X_norm <- clr_rows(X_filt, pseudocount = PSEUDOCOUNT)
} else if (USE_METHOD == "CSS") {
  msg("[4] Normalization = CSS (metagenomeSeq)")
  X_norm <- css_rows(X_filt)
  # optional log1p
  X_norm <- log1p(X_norm)
} else {
  stop("USE_METHOD must be 'CLR' or 'CSS'")
}

# scale features for PCA
X_scaled <- scale(X_norm)

## -----------------------------
## 5) PCA (linear)
## -----------------------------
N_PCS = 20

pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
npc <- min(N_PCS, ncol(pca$x))
PC <- pca$x[, 1:npc, drop = FALSE]
rownames(PC) <- rownames(X_scaled)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pve <- data.frame(PC = paste0("PC", seq_along(var_exp)),
                  variance_explained = var_exp)

fwrite(pve, file.path(OUTDIR, "Tables", "PCA_variance_explained.tsv"), sep="\t")
fwrite(data.frame(MAG=rownames(PC), PC, check.names = FALSE),
       file.path(OUTDIR, "Tables", "PCA_scores.tsv"), sep="\t")

## -----------------------------
## 6) UMAP / tSNE on PCs
## -----------------------------
um <- uwot::umap(PC, n_neighbors = UMAP_NN, min_dist = UMAP_MD, metric = "euclidean")
colnames(um) <- c("UMAP1","UMAP2")
rownames(um) <- rownames(PC)

ts <- Rtsne::Rtsne(PC, perplexity = min(TSNE_PERP, floor((nrow(PC)-1)/3)),
                   check_duplicates = FALSE)$Y
colnames(ts) <- c("TSNE1","TSNE2")
rownames(ts) <- rownames(PC)

## -----------------------------
## 7) Distance for evaluation / PCoA
## Bray-Curtis on corrected counts (X_corr or X_filt)
## -----------------------------
bc <- vegan::vegdist(X_filt, method = "bray")
# PCoA
pcoa <- cmdscale(bc, k = 2, eig = TRUE)
pcoa_xy <- pcoa$points
colnames(pcoa_xy) <- c("PCoA1","PCoA2")
rownames(pcoa_xy) <- rownames(X_filt)

## -----------------------------
## 8) Clustering (no preset K): HDBSCAN on PC space
## -----------------------------

minpts_grid <- c(5, 10, 15, 20, 30, 40, 50)

safe_mean_prob <- function(h){
  p <- h$probabilities
  if (is.null(p)) return(NA_real_)
  if (!is.numeric(p)) p <- suppressWarnings(as.numeric(p))
  if (all(is.na(p))) return(NA_real_)
  sel <- h$cluster != 0
  if (!any(sel)) return(NA_real_)
  mean(p[sel], na.rm = TRUE)
}

safe_mean_persistence <- function(h){
  # persistence 通常是每个 cluster 的稳定性分数
  pers <- h$cluster_scores
  if (is.null(pers)) return(NA_real_)
  if (!is.numeric(pers)) pers <- suppressWarnings(as.numeric(pers))
  if (all(is.na(pers))) return(NA_real_)
  mean(pers, na.rm = TRUE)
}

res <- bind_rows(lapply(minpts_grid, function(m){
  h <- dbscan::hdbscan(PC, minPts = m)
  
  data.frame(
    minPts = m,
    n_cluster  = length(setdiff(unique(h$cluster), 0)),
    noise_frac = mean(h$cluster == 0),
    mean_prob  = safe_mean_prob(h),
    mean_persistence = safe_mean_persistence(h)
  )
}))

res



HDBSCAN_MINPTS = 20
msg("[8] HDBSCAN on PC space (minPts=%d)", HDBSCAN_MINPTS)
hdb <- dbscan::hdbscan(PC, minPts = HDBSCAN_MINPTS)
cl_hdb <- hdb$cluster
# 0 = noise
n_clu <- length(setdiff(unique(cl_hdb), 0))
noise_frac <- mean(cl_hdb == 0)

msg("    HDBSCAN clusters=%d, noise=%.2f", n_clu, noise_frac)

# Optional: DBSCAN (eps chosen from kNN distance knee can be complicated; here just export kNNdist)
if (RUN_DBSCAN) {
  k_eps <- 10
  knn_d <- dbscan::kNNdist(PC, k = k_eps)
  # save distances for you to choose eps (look for elbow)
  fwrite(data.frame(knn_dist=sort(knn_d)),
         file.path(OUTDIR, "Tables", paste0("DBSCAN_kNNdist_k",k_eps,".tsv")),
         sep="\t")
  msg("    [DBSCAN] wrote kNNdist table (choose eps by elbow): Tables/DBSCAN_kNNdist_k%d.tsv", k_eps)
}

## -----------------------------
## 9) Quality metrics (no K selection, just QC)
## -----------------------------
# silhouette on PC euclidean distance, excluding noise
PC_dist <- dist(PC)
sil_hdb <- NA_real_
if (n_clu >= 2) {
  keep_idx <- cl_hdb != 0
  if (sum(keep_idx) > 10 && length(unique(cl_hdb[keep_idx])) >= 2) {
    sil_hdb <- silhouette_score(cl_hdb[keep_idx], dist(PC[keep_idx, , drop=FALSE]))
  }
}
msg("[9] silhouette (HDBSCAN, excluding noise) = %s", ifelse(is.na(sil_hdb), "NA", sprintf("%.3f", sil_hdb)))

qc <- data.frame(
  method = "HDBSCAN",
  minPts = HDBSCAN_MINPTS,
  n_clusters = n_clu,
  noise_frac = noise_frac,
  silhouette_excl_noise = sil_hdb
)
fwrite(qc, file.path(OUTDIR, "Tables", "clustering_QC.tsv"), sep="\t")

## -----------------------------
## 10) Output tables
## -----------------------------
res_tbl <- data.frame(
  MAG = rownames(PC),
  comp_frac = comp[rownames(PC)],
  hdbscan_cluster = cl_hdb,
  UMAP1 = um[,1], UMAP2 = um[,2],
  TSNE1 = ts[,1], TSNE2 = ts[,2],
  PCoA1 = pcoa_xy[rownames(PC), 1],
  PCoA2 = pcoa_xy[rownames(PC), 2],
  stringsAsFactors = FALSE
)

fwrite(res_tbl, file.path(OUTDIR, "Tables", "MAG_embeddings_clusters.tsv"), sep="\t")

## -----------------------------
## 11) Plots
## -----------------------------
theme_sc <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none',
        axis.ticks.length = unit(2,'mm'))

# helper for discrete color with noise
res_tbl <- res_tbl %>%
  mutate(hdbscan_cluster = ifelse(hdbscan_cluster == 0, "Noise", as.character(hdbscan_cluster)),
         hdbscan_cluster = factor(hdbscan_cluster, levels = c("Noise", sort(setdiff(unique(as.character(hdb$cluster)), "0")))))


library(Polychrome)

n_cluster <- nlevels(res_tbl$hdbscan_cluster) - 1  # 不含 noise

# 为 cluster 生成高区分度颜色
set.seed(123)
cluster_cols <- Polychrome::createPalette(
  n_cluster,
  c("#3B4CC0", "#B40426"),  # 冷→暖，避免荧光
  M = 5000
)

# 最终调色表
cluster_palette <- c(
  "Noise" = "grey80",
  setNames(cluster_cols, seq_len(n_cluster))
)


p_umap <- ggplot(res_tbl, aes(UMAP1, UMAP2, color = hdbscan_cluster)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = cluster_palette) +
  theme_sc +
  labs(
    title = sprintf("UMAP (PCA→UMAP) + HDBSCAN (minPts=%d)", HDBSCAN_MINPTS),
    color = "Cluster"
  )
p_umap

p_tsne <- ggplot(res_tbl, aes(TSNE1, TSNE2, color = hdbscan_cluster)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = cluster_palette) +
  theme_sc
  
p_tsne

p_pcoa <- ggplot(res_tbl, aes(PCoA1, PCoA2, color = hdbscan_cluster)) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = cluster_palette) +
  theme_sc +
  labs(
    title = "PCoA (Bray–Curtis on corrected copy numbers) + HDBSCAN labels",
    color = "Cluster"
  )
p_pcoa

p_pca <- ggplot(
  data.frame(
    PC1 = PC[,1],
    PC2 = PC[,2],
    cluster = res_tbl$hdbscan_cluster
  ),
  aes(PC1, PC2, color = cluster)
) +
  geom_point(size = 1.2, alpha = 0.85) +
  scale_color_manual(values = cluster_palette) +
  theme_sc +
  labs(
    title = "PCA (CLR/CSS) + HDBSCAN labels",
    color = "Cluster"
  )
p_pca


ggsave(file.path(OUTDIR, "Figures", "UMAP_HDBSCAN.pdf"), p_umap, width = 5.2, height = 5.2)
ggsave(file.path(OUTDIR, "Figures", "TSNE_HDBSCAN-1.pdf"), p_tsne, width = 6, height = 6)
ggsave(file.path(OUTDIR, "Figures", "PCoA_Bray_HDBSCAN.pdf"), p_pcoa, width = 5.2, height = 5.2)
ggsave(file.path(OUTDIR, "Figures", "PCA_HDBSCAN.pdf"), p_pca, width = 5.2, height = 5.2)

msg("DONE.
Outputs:
- Tables/MAGxFeature_copy_corrected.tsv
- Tables/feature_filtering_summary.tsv
- Tables/PCA_scores.tsv
- Tables/MAG_embeddings_clusters.tsv
- Tables/clustering_QC.tsv
- Figures/UMAP_HDBSCAN.pdf
- Figures/TSNE_HDBSCAN.pdf
- Figures/PCoA_Bray_HDBSCAN.pdf
- Figures/PCA_HDBSCAN.pdf
")

#################
#MAG 大类功能鉴定
###################

#################
# MAG cluster functional definition by KEGG Modules
#################

msg("[12] Cluster functional definition using KEGG Modules ...")

# ---- config (module) ----
IN_MODULE <- "./微生物/MAGxKEGG_Module_copy_number.tsv"   # 你拿到的矩阵
TOPN_MODULES <- 25      # 每个 cluster 输出多少个显著 module
MIN_PREV_IN_CLUSTER <- 0.25  # module 在 cluster 内至少 25% MAG 存在才考虑（避免单个MAG驱动）
MIN_ABS_LOG2FC <- 2        # 差异阈值（可调）
ALPHA_FDR <- 0.05

# ---- 1) read module matrix ----
mod <- fread(IN_MODULE, data.table = FALSE)
colnames(mod)[1] <- "MAG"
mod$MAG <- strip_ext(as.character(mod$MAG))

# ---- 1b) completeness correction for Module copy numbers ----
# 说明：你的 comp_frac 是 0-1 之间的比例（已在前面 meta 处理过）
# 校正：copy_corrected = copy / comp_frac

# 先把 comp_frac 合并进来（只保留聚类用到的 MAG）
mod <- mod %>%
  inner_join(res_tbl %>% select(MAG, comp_frac, hdbscan_cluster), by = "MAG")

mod_cols <- setdiff(colnames(mod), c("MAG", "comp_frac", "hdbscan_cluster"))

# 数值化
mod <- mod %>% mutate(across(all_of(mod_cols), ~ suppressWarnings(as.numeric(.x))))
mod[is.na(mod)] <- 0

# 校正（按行除以 comp_frac）
comp_vec <- mod$comp_frac
comp_vec[!is.finite(comp_vec) | comp_vec <= 0] <- NA_real_

M_raw <- as.matrix(mod[, mod_cols, drop = FALSE])
M_corr <- sweep(M_raw, 1, comp_vec, `/`)
M_corr[!is.finite(M_corr)] <- 0

# 用校正后的值替换回 mod（后面差异分析就自动用 corrected）
mod[, mod_cols] <- M_corr

# 可选：把校正后的模块矩阵保存，便于复用/审稿复现
fwrite(
  data.frame(MAG = mod$MAG, M_corr, check.names = FALSE),
  file.path(OUTDIR, "Tables", "MAGxKEGG_Module_copy_corrected.tsv"),
  sep = "\t"
)

msg("[12] Module copy numbers corrected by completeness (copy/comp_frac).")


# 把 cluster 名称标准化：Noise / C1..C18（你前面 res_tbl 已经 factor 了）
# 这里保证字符一致
mod$cluster <- as.character(mod$hdbscan_cluster)

# ---- 2) helper: per-cluster stats & enrichment test ----
# 我们对每个 module 做 “cluster vs others” 的 Wilcoxon（稳健、非正态友好）
# 并输出：
#   prev_in_cluster / prev_in_others
#   mean_in_cluster / mean_in_others
#   log2FC = log2((mean_in_cluster+eps)/(mean_in_others+eps))
#   pvalue, FDR

calc_cluster_enrich <- function(df_long, clu, eps = 1e-6){
  d1 <- df_long %>% filter(cluster == clu)
  d0 <- df_long %>% filter(cluster != clu)  # Others：排除 Noise（你也可把 Noise 算进 Others）
  if (nrow(d1) < 5 || nrow(d0) < 20) return(NULL)
  
  res <- d1 %>%
    group_by(Module) %>%
    summarise(
      mean_in = mean(copy, na.rm = TRUE),
      prev_in = mean(copy > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      d0 %>%
        group_by(Module) %>%
        summarise(
          mean_out = mean(copy, na.rm = TRUE),
          prev_out = mean(copy > 0, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Module"
    ) %>%
    mutate(
      log2FC = log2((mean_in + eps) / (mean_out + eps))
    )
  
  # Wilcoxon pvalue（逐 module）
  # 注意：对 1497 MAG × 上千 module，会比较慢；但能跑。要更快可以改成稀疏矩阵 + fast test
  pvals <- sapply(res$Module, function(m){
    x <- d1$copy[d1$Module == m]
    y <- d0$copy[d0$Module == m]
    # 若两边全 0 或变化太少，给 NA
    if (all(x == 0) && all(y == 0)) return(NA_real_)
    # wilcox 至少要有变异
    if (length(unique(c(x, y))) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(x, y)$p.value)
  })
  
  res$pvalue <- pvals
  res$fdr <- p.adjust(res$pvalue, method = "BH")
  res$cluster <- clu
  res
}

# ---- 3) convert to long table ----
mod_long <- mod %>%
  select(MAG, cluster, all_of(mod_cols)) %>%
  pivot_longer(cols = all_of(mod_cols),
               names_to = "Module",
               values_to = "copy") %>%
  mutate(
    copy = as.numeric(copy),
    copy = ifelse(is.na(copy), 0, copy)
  )

# 基础统计表：cluster 内 module prevalence / mean（不做差异）
cluster_module_stats <- mod_long %>%
  group_by(cluster, Module) %>%
  summarise(
    n_MAG = n(),
    prev = mean(copy > 0, na.rm = TRUE),
    mean_copy = mean(copy, na.rm = TRUE),
    median_copy = median(copy, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(cluster_module_stats,
       file.path(OUTDIR, "Tables", "Cluster_Module_stats.tsv"),
       sep = "\t")

# ---- 4) enrichment per cluster (exclude Noise) ----
clusters_use <- unique(mod$cluster)

enrich_list <- lapply(clusters_use, function(clu){
  msg("    [enrich] cluster=%s", clu)
  calc_cluster_enrich(mod_long, clu)
})
enrich <- bind_rows(enrich_list)

# 过滤：cluster 内至少一定 prevalence，且差异达到阈值
enrich_filt <- enrich %>%
  filter(!is.na(fdr)) %>%
  filter(prev_in >= MIN_PREV_IN_CLUSTER) %>%
  filter(abs(log2FC) >= MIN_ABS_LOG2FC) %>%
  arrange(cluster, fdr)


enrich_filt <- left_join(enrich_filt, module_name, by='Module')

fwrite(enrich,
       file.path(OUTDIR, "Tables", "Cluster_Module_enrichment_all.tsv"),
       sep = "\t")
fwrite(enrich_filt,
       file.path(OUTDIR, "Tables", "Cluster_Module_enrichment.tsv"),
       sep = "\t")

# ---- 5) auto-generate functional definition (draft) ----
module_name <- fread('./微生物/MAG_COG_cluster/Tables/module name.txt') %>%
  distinct(Module, name)


# 1) 确保 Module 是 character，且 module_name 也去重
module_name2 <- module_name %>%
  mutate(Module = as.character(Module)) %>%
  distinct(Module, name)

enrich_filt2 <- enrich_filt %>%
  mutate(Module = as.character(Module)) %>%     # 关键：保证 join 键一致
  left_join(module_name2, by = "Module", suffix = c("", ".kegg"))

# 2) 兼容：join 后名字列可能叫 name / name.kegg / name.x / name.y
name_candidates <- intersect(c("name", "name.kegg", "name.x", "name.y"), colnames(enrich_filt2))
if (length(name_candidates) == 0) {
  stop("Join 后仍找不到模块名称列（name/name.kegg/name.x/name.y）。请检查 enrich_filt 的 Module 列名是否正确。")
}

# 3) 统一模块名称列到 module_name（优先用 join 带来的 name.kegg，其次用原来的 name）
#    注意：如果 enrich_filt 本来就有 name（比如你截图里 enrich 有 name），那 join 的会是 name.kegg
enrich_filt2 <- enrich_filt2 %>%
  mutate(
    module_name = coalesce(!!!rlang::syms(name_candidates)),  # 取第一个非 NA 的
    module_name = ifelse(is.na(module_name) | module_name == "", Module, module_name),
    mod_id_name = paste0(Module, ": ", module_name)
  )

collapse_unique <- function(x, sep=";"){
  x <- x[!is.na(x) & x != ""]
  paste(unique(x), collapse = sep)
}

top_mods <- enrich_filt2 %>%
  filter(fdr <= ALPHA_FDR, log2FC > 2) %>%
  group_by(cluster) %>%
  arrange(fdr, desc(log2FC)) %>%
  slice_head(n = TOPN_MODULES) %>%
  summarise(
    n_sig = n(),
    top_modules_id = collapse_unique(Module),
    top_modules_name = collapse_unique(module_name),
    top_modules_id_name = collapse_unique(mod_id_name),
    .groups = "drop"
  )

fallback_mods <- enrich_filt2 %>%
  group_by(cluster) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = TOPN_MODULES) %>%
  summarise(
    top_modules_id_fallback = collapse_unique(Module),
    top_modules_name_fallback = collapse_unique(module_name),
    top_modules_id_name_fallback = collapse_unique(mod_id_name),
    .groups = "drop"
  )

cluster_def <- data.frame(cluster = clusters_use, stringsAsFactors = FALSE) %>%
  left_join(top_mods, by = "cluster") %>%
  left_join(fallback_mods, by = "cluster") %>%
  mutate(
    definition_draft = ifelse(!is.na(top_modules_id_name) & n_sig > 0,
                              paste0("Enriched KEGG modules: ", top_modules_id_name),
                              paste0("Top modules (no FDR<", ALPHA_FDR, "): ", top_modules_id_name_fallback))
  )

fwrite(cluster_def,
       file.path(OUTDIR, "Tables", "Cluster_Function_Definition.tsv"),
       sep = "\t")





