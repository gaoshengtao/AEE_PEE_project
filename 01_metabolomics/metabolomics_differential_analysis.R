############################################################
# AEE vs PEE non-targeted metabolomics differential analysis
# Method: OPLS-DA (1+1) + Wilcoxon p<0.05
# Input:
#   1) 归一化代谢组数据.xlsx
#   2) Supplemental file 1-VFA_rumenotype_assignments.xlsx
# Output:
#   - results_all_metabolites.tsv
#   - results_sig_OPLSDA1_Wilcox005.tsv
#   - OPLSDA_scoreplot.pdf / VIP_top20.pdf / Volcano.pdf
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# ---------- 0) paths ----------
assign_xlsx  <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx"
metabo_xlsx  <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/代谢组数据/归一化代谢组数据.xlsx"
assign_sheet <- "VFA_rumenotype_assignments"

outdir <- "AEE_PEE_metabolomics_OPLSDA_Wilcox"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1) read data ----------
met0 <- readxl::read_excel(metabo_xlsx, sheet = 1) %>% as.data.frame()
ass0 <- readxl::read_excel(assign_xlsx, sheet = assign_sheet) %>% as.data.frame()

# keep only needed columns in assignment
ass <- ass0 %>%
  transmute(
    SampleID   = as.character(SampleID),
    Rumenotype = as.character(Rumenotype)
  ) %>%
  filter(!is.na(SampleID), !is.na(Rumenotype))

# map to AEE/PEE labels (Acetate_type -> AEE; Propionate_type -> PEE)
ass <- ass %>%
  mutate(
    Group = dplyr::case_when(
      Rumenotype == "Acetate_type"     ~ "AEE",
      Rumenotype == "Propionate_type"  ~ "PEE",
      TRUE                             ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# ---------- 2) identify sample columns in metabolomics table ----------
# sample columns look like Dxxxx / Dxxxxx; we intersect with assignment SampleID
sample_cols <- intersect(colnames(met0), ass$SampleID)
if (length(sample_cols) < 6) {
  stop("Matched sample columns too few. Check SampleID naming consistency between two files.")
}

# annotation columns = non-sample columns
anno_cols <- setdiff(colnames(met0), sample_cols)

anno_cols <- anno_cols[1:25]
# ---------- 3) build feature table ----------
# create a robust feature_id for each metabolite row
met <- met0 %>%
  mutate(
    feature_id = if ("ID" %in% colnames(met0)) as.character(ID) else paste0("feat_", row_number()),
    feature_id = make.unique(feature_id),
    Metabolite = if ("Metabolite" %in% colnames(met0)) as.character(Metabolite) else NA_character_,
    mz         = if ("m/z" %in% colnames(met0)) as.numeric(`m/z`) else NA_real_,
    rt         = if ("Retention.time" %in% colnames(met0)) as.numeric(Retention.time) else NA_real_,
    Mode       = if ("Mode" %in% colnames(met0)) as.character(Mode) else NA_character_
  )

# matrix: samples x features
X_feat_by_samp <- met %>%
  select(feature_id, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric))

# handle missing values: median impute per feature (row)
mat <- as.matrix(X_feat_by_samp[, sample_cols, drop = FALSE])
rownames(mat) <- X_feat_by_samp$feature_id

row_median <- function(x) median(x, na.rm = TRUE)
meds <- apply(mat, 1, row_median)
for (i in seq_len(nrow(mat))) {
  nas <- is.na(mat[i, ])
  if (any(nas)) mat[i, nas] <- meds[i]
}

# transpose to ropls format: samples x variables
X <- t(mat)
# align Y
meta <- ass %>%
  filter(SampleID %in% rownames(X)) %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  arrange(match(SampleID, rownames(X)))

X <- X[meta$SampleID, , drop = FALSE]
Y <- factor(meta$Group, levels = c("AEE", "PEE"))

# optional: drop near-zero variance features
nzv <- apply(X, 2, sd, na.rm = TRUE)
X <- X[, nzv > 1e-8, drop = FALSE]

# ---------- 4) OPLS-DA (ropls) ----------
# Install if missing:
if (!requireNamespace("ropls", quietly = TRUE)) {
  message("Package 'ropls' not found. Installing from Bioconductor...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ropls")
}
suppressPackageStartupMessages(library(ropls))

set.seed(123)

# OPLS-DA 1+1 (1 predictive + 1 orthogonal)
# permI 设为 0/200 都行；你先跑 200 看稳不稳
opls_fit <- opls(
  x = X,
  y = Y,
  predI = NA,
  orthoI = NA,
  permI = 200
)


# extract scores
score_df <- data.frame(
  SampleID = rownames(X),
  Group = Y,
  t1 = opls_fit@scoreMN[, 1],
  to1 = if (!is.null(opls_fit@orthoScoreMN) && ncol(opls_fit@orthoScoreMN) >= 1) opls_fit@orthoScoreMN[, 1] else 0
)

# extract VIP
vip <- NULL
if (!is.null(opls_fit@vipVn)) {
  vip <- opls_fit@vipVn
} else if ("getVipVn" %in% getNamespaceExports("ropls")) {
  vip <- ropls::getVipVn(opls_fit)
} else {
  warning("VIP not found in ropls object. Will proceed without VIP filtering.")
}
vip_df <- data.frame(feature_id = names(vip), VIP = as.numeric(vip), row.names = NULL)

# ---------- 5) Wilcoxon per feature ----------
# compute median difference (AEE - PEE) on your normalized scale
XA <- X[Y == "AEE", , drop = FALSE]
XP <- X[Y == "PEE", , drop = FALSE]

medA <- apply(XA, 2, mean, na.rm = TRUE)
medP <- apply(XP, 2, mean, na.rm = TRUE)
delta_med <- medA - medP                # 类似 logFC（如果你的数据本身就是log/CLR）
log2fc <- log2(medA / medP)
fc_like   <- 2^(delta_med)              # 仅作为直观 fold change（若数据并非log2则仅参考）

wilcox_p <- sapply(colnames(X), function(feat) {
  suppressWarnings(wilcox.test(XA[, feat], XP[, feat])$p.value)
})

res <- data.frame(
  feature_id = colnames(X),
  median_AEE = as.numeric(medA[colnames(X)]),
  median_PEE = as.numeric(medP[colnames(X)]),
  delta_med  = as.numeric(delta_med[colnames(X)]),
  log2fc = as.numeric(log2fc[colnames(X)]),
  FC_like    = as.numeric(fc_like[colnames(X)]),
  p_wilcox   = as.numeric(wilcox_p[colnames(X)]),
  stringsAsFactors = FALSE
) %>%
  mutate(padj = p.adjust(p_wilcox, method = "BH"))

# merge annotation + VIP
anno_df <- met %>%
  select(all_of(c("feature_id", intersect(anno_cols, colnames(met))))) %>%
  distinct(feature_id, .keep_all = TRUE)

res_all <- res %>%
  left_join(vip_df, by = "feature_id") %>%
  left_join(anno_df, by = "feature_id")

# ---------- 6) significance filter ----------
use_vip_filter <- TRUE      # 你如果只要 “oplsda 1 + wilcox 0.05” 不想加VIP，就改 FALSE
vip_cut <- 1.0
p_cut <- 0.05

res_sig <- res_all %>%
  filter(p_wilcox < p_cut) %>%
  { if (use_vip_filter && "VIP" %in% colnames(.)) filter(., VIP >= vip_cut) else . } %>%
  arrange(p_wilcox)

res_sig$Dir <- ifelse(res_sig$median_AEE > res_sig$median_PEE, 'AEE', 'PEE')

# ---------- 7) output tables ----------
readr::write_tsv(res_all, file.path(outdir, "results_all_metabolites.tsv"))
readr::write_tsv(res_sig, file.path(outdir, "results_sig_OPLSDA1_Wilcox005.tsv"))

# ---------- 8) plots ----------
# (1) OPLS-DA score plot: t[1] vs to[1]
vol_color <- c('AEE' = '#2e528f', 'PEE' = '#b41d23', 'N.S.' = '#efefef')

p_score <- ggplot(score_df, aes(x = t1, y = to1, color = Group)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    x = "t[1] (predictive)",
    y = "to[1] (orthogonal)"
  )+
  theme(panel.grid = element_blank(),
        # axis.text = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.9),
        axis.ticks.length = unit(2,'mm')
        )+
  scale_color_manual(values = vol_color)
p_score
ggsave(file.path(outdir, "OPLSDA_scoreplot.pdf"), p_score, width = 5, height = 5)


# ------------------------------------------------------------
# 1) 提取 R2X / R2Y / Q2
# ------------------------------------------------------------
perf_df <- opls_fit@summaryDF

# ropls 的 summaryDF 有时行名是 "Total"，有时是 "overall"/"Model"（不同版本可能略有差异）
# 这里做个稳健提取：
get_metric <- function(df, row_candidates = c("Total","overall","Model","All"), colname){
  rn <- rownames(df)
  row_pick <- intersect(row_candidates, rn)
  if (length(row_pick) == 0) {
    stop("Cannot find Total/overall row in opls_fit@summaryDF. Please check structure(perf_df).")
  }
  as.numeric(df[row_pick[1], colname])
}

R2X <- get_metric(perf_df, colname = "R2X(cum)")
R2Y <- get_metric(perf_df, colname = "R2Y(cum)")
Q2  <- get_metric(perf_df, colname = "Q2(cum)")

cat("R2X =", R2X, " R2Y =", R2Y, " Q2 =", Q2, "\n")

perf_out <- tibble(
  R2X_cum = R2X,
  R2Y_cum = R2Y,
  Q2_cum  = Q2
)
write_tsv(perf_out, file.path(outdir, "OPLS_performance_R2X_R2Y_Q2.tsv"))

# ------------------------------------------------------------
# 2) Permutation test 图 + 原始矩阵导出
# ------------------------------------------------------------
pdf(file.path(outdir, "OPLS_perm_plot.pdf"), width = 6, height = 5)
plot(opls_fit, typeVc = "permutation")
dev.off()

# permutation 原始矩阵（用于附表/自己检查）
perm_mat <- NULL
if (!is.null(opls_fit@suppLs) && "permMN" %in% names(opls_fit@suppLs)) {
  perm_mat <- opls_fit@suppLs[["permMN"]]
  perm_df <- as.data.frame(perm_mat) %>% tibble::rownames_to_column("perm_id")
  write_tsv(perm_df, file.path(outdir, "OPLS_perm_matrix.tsv"))
}

# ------------------------------------------------------------
# 3) CV-ANOVA p-value（模型显著性检验）
# ------------------------------------------------------------
cv_anova_p <- NA
if (!is.null(opls_fit@suppLs) && "cv-anova" %in% names(opls_fit@suppLs)) {
  cv_anova_p <- opls_fit@suppLs[["cv-anova"]]
}
cat("CV-ANOVA p =", cv_anova_p, "\n")

write_tsv(
  tibble(cv_anova_p = as.numeric(cv_anova_p)),
  file.path(outdir, "OPLS_CV_ANOVA_p.tsv")
)

# ------------------------------------------------------------
# 4) x-score plot（score plot）—— ropls 自带
#    注意：你之前 ggplot 的 t1/to1 也可以；这个是 ropls 官方默认样式
# ------------------------------------------------------------
pdf(file.path(outdir, "OPLS_score_plot_ropls_xscore.pdf"), width = 6, height = 5)
plot(opls_fit, typeVc = "x-score")
dev.off()

# ------------------------------------------------------------
# 5) VIP 图（Top 20）—— 你之前已有，这里再给一个 ropls 提取 + 画法
# ------------------------------------------------------------
vip <- NULL
if (!is.null(opls_fit@vipVn)) {
  vip <- opls_fit@vipVn
} else if ("getVipVn" %in% getNamespaceExports("ropls")) {
  vip <- ropls::getVipVn(opls_fit)
}

if (!is.null(vip)) {
  vip_df <- tibble(feature_id = names(vip), VIP = as.numeric(vip)) %>%
    arrange(desc(VIP))
  
  write_tsv(vip_df, file.path(outdir, "OPLS_VIP_all.tsv"))
  
  vip_top20 <- vip_df %>% slice_head(n = 20) %>%
    mutate(feature_id = factor(feature_id, levels = rev(feature_id)))
  
  p_vip <- ggplot(vip_top20, aes(x = feature_id, y = VIP)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(title = "Top 20 VIP (OPLS-DA)", x = NULL, y = "VIP")
  
  ggsave(file.path(outdir, "OPLS_VIP_top20.pdf"), p_vip, width = 8, height = 6)
}

# ------------------------------------------------------------
# 6) S-plot：p vs p(corr)
#    p = loading（与 predictive component 方向相关的权重）
#    p(corr) = 变量与 t1 score 的相关系数（反映稳定性/相关方向）
# ------------------------------------------------------------

# 6.1 数据标准化（S-plot 建议用 scaled 数据）
# ropls内部会做 scaling，但为了你可控、可复现，这里手动 scale 一份
X_scaled <- scale(X)

# 6.2 提取 loadings（p）
# getLoadingMN 返回 loadings matrix：变量 × 成分
loadingMN <- ropls::getLoadingMN(opls_fit)
loading_p <- loadingMN[, 1]    # 第1个 predictive component 的 loading（p）

# 6.3 提取 scores（t1）
scoreMN <- ropls::getScoreMN(opls_fit)
scores_t1 <- scoreMN[, 1]

# 6.4 计算 p(corr)：每个变量与 t1 的相关系数
pcorr <- apply(X_scaled, 2, function(v) suppressWarnings(cor(v, scores_t1)))

splot_df <- tibble(
  feature_id = names(loading_p),
  p = as.numeric(loading_p),
  pcorr = as.numeric(pcorr)
)

splot_df <- left_join(splot_df, res_all,by='feature_id')

splot_df$color <- ifelse((splot_df$median_AEE > splot_df$median_PEE) & splot_df$VIP > 1 & splot_df$p_wilcox, 'AEE', 
                         ifelse((splot_df$median_AEE < splot_df$median_PEE) & splot_df$VIP > 1 & splot_df$p_wilcox, 'PEE', 'N.S.')
                         )

write_tsv(splot_df, file.path(outdir, "OPLS_Splot_data.tsv"))

# 6.5 画 S-plot
p_splot <- ggplot(splot_df, aes(x = p, y = pcorr, color = color)) +
  geom_point(size = 1.2, alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "S-plot (p vs p(corr))", x = "p (loading)", y = "p(corr) (corr with t[1])") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = vol_color)


p_splot

ggsave(file.path(outdir, "OPLS_Splot.pdf"), p_splot, width = 6.5, height = 5.5)

# ------------------------------------------------------------
# 7) Loading plot（可选，展示 p 值最大的变量）
# ------------------------------------------------------------
topN <- 30
loading_top <- splot_df %>%
  mutate(abs_p = abs(p)) %>%
  arrange(desc(abs_p)) %>%
  slice_head(n = topN) %>%
  mutate(feature_id = factor(feature_id, levels = rev(feature_id)))

p_loading <- ggplot(loading_top, aes(x = feature_id, y = p)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = paste0("Top ", topN, " loadings (p) on t[1]"), x = NULL, y = "p (loading)")

ggsave(file.path(outdir, paste0("OPLS_Loading_top", topN, ".pdf")), p_loading, width = 8, height = 7)

cat("\n[Done] Extra OPLS-DA figures saved in: ", normalizePath(outdir), "\n")


# (3) Volcano: delta_med vs -log10(p)
res_plot <- res_all %>%
  mutate(
    neglog10p = -log10(p_wilcox),
    Sig = ifelse(p_wilcox < p_cut & VIP >1 , "p<0.05", "NS"),
    Dir = ifelse(p_wilcox < p_cut & VIP >1 & median_AEE > median_PEE, 'AEE',
                 ifelse(p_wilcox < p_cut & VIP >1 & median_AEE < median_PEE, 'PEE','N.S.'
                 ))
  )

vol_color <- c('AEE'='#2e528f', 'PEE'='#b41d23', 'N.S.'= '#efefef')

p_vol <- ggplot(res_plot, aes(x = log2fc, y = -log10(p_wilcox), colour = Dir)) +
  geom_point( size = 2) +
  theme_bw(base_size = 12) +
  # labs( x = "delta median (AEE - PEE)", y = "-log10(Wilcoxon p)") +
  geom_hline(yintercept = -log10(p_cut), linetype = 2) +
  scale_color_manual(values = vol_color)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        axis.ticks.length = unit(2,'mm'))

p_vol

table(res_plot$Dir)

ggsave(file.path(outdir, "Volcano.pdf"), p_vol, width = 5, height = 5)

# ---------- 9) quick summary ----------
cat("\n================ SUMMARY ================\n")
cat("Total matched samples:", nrow(X), "\n")
cat("Total features used:", ncol(X), "\n")
cat("Wilcoxon p<0.05 features:", sum(res_all$p_wilcox < 0.05, na.rm = TRUE), "\n")
if ("VIP" %in% colnames(res_all)) {
  cat("With VIP>=", vip_cut, "and p<0.05 features:", nrow(res_sig), "\n")
} else {
  cat("VIP unavailable; significant features only by p<0.05:", nrow(res_sig), "\n")
}
cat("Outputs saved in:", normalizePath(outdir), "\n")
cat("========================================\n")



############################################################
# AEE vs PEE non-targeted metabolomics
# PCoA + adonis2 (PERMANOVA)
# Input:
#   1) 归一化代谢组数据.xlsx
#   2) Supplemental file 1-VFA_rumenotype_assignments.xlsx (sheet: VFA_rumenotype_assignments)
# Output:
#   - PCoA_plot.pdf
#   - adonis2_result.tsv
#   - distance_matrix.rds
############################################################
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(vegan)    # vegdist + adonis2
  library(ggplot2)
})


# ---------- 4) distance + PCoA ----------
# Choose distance method:
# - "bray" for abundance/intensity-like data (most common)
# - "euclidean" if already log/CLR and you want Euclidean geometry
distance_method <- "bray"

# vegdist requires non-negative for bray; if your data has negatives (CLR), switch to euclidean
if (distance_method == "bray" && any(X < 0, na.rm = TRUE)) {
  warning("Detected negative values; Bray-Curtis is invalid. Switching to Euclidean.")
  distance_method <- "euclidean"
}

dist_obj <- vegan::vegdist(X, method = distance_method)
saveRDS(dist_obj, file.path(outdir, "distance_matrix.rds"))

# PCoA via cmdscale
pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE, add = TRUE)
pcoa_df <- data.frame(
  SampleID = rownames(X),
  Group    = meta$Group,
  PC1      = pcoa$points[, 1],
  PC2      = pcoa$points[, 2]
)

var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
pc1_pct <- round(100 * var_expl[1], 2)
pc2_pct <- round(100 * var_expl[2], 2)

# ---------- 5) adonis2 (PERMANOVA) ----------
set.seed(123)
adon <- vegan::adonis2(dist_obj ~ Group, data = meta, permutations = 999, by = "terms")
adon_tbl <- as.data.frame(adon) %>%
  tibble::rownames_to_column("Term")

readr::write_tsv(adon_tbl, file.path(outdir, "adonis2_result.tsv"))

# also compute beta-dispersion (optional but recommended)
# If dispersion differs a lot, PERMANOVA may reflect dispersion not centroid
bd <- vegan::betadisper(dist_obj, meta$Group)
bd_anova <- anova(bd)
bd_perm  <- permutest(bd, permutations = 999)

bd_out <- list(
  betadisper_anova = bd_anova,
  betadisper_permutest = bd_perm
)
saveRDS(bd_out, file.path(outdir, "betadisper_check.rds"))

# ---------- 6) plot ----------
# Build annotation text from adonis2
pval <- adon_tbl$`Pr(>F)`[adon_tbl$Term == "Group"]
R2   <- adon_tbl$R2[adon_tbl$Term == "Group"]

lab_txt <- paste0(
  "adonis2: R2 = ", sprintf("%.3f", R2),
  ", p = ", format.pval(pval, digits = 3, eps = 1e-3),
  "\nDistance: ", distance_method
)

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCoA of non-targeted metabolomics (AEE vs PEE)",
    x = paste0("PCoA1 (", pc1_pct, "%)"),
    y = paste0("PCoA2 (", pc2_pct, "%)")
  ) 
  
  p

ggsave(file.path(outdir, "PCoA_plot.pdf"), p, width = 6.8, height = 5.6)

# ---------- 7) summary ----------
cat("\n================ SUMMARY ================\n")
cat("Samples:", nrow(X), "\n")
cat("Features:", ncol(X), "\n")
cat("Distance:", distance_method, "\n")
cat("adonis2 R2:", sprintf("%.4f", R2), "\n")
cat("adonis2 p:", format.pval(pval, digits = 3, eps = 1e-3), "\n")
cat("Outputs:", normalizePath(outdir), "\n")
cat("========================================\n")



res_sig$dir <- ifelse(res_sig$median_AEE - res_sig$median_PEE, 'AEE', 'PEE')

meta <- arrange(meta, Group)
res_sig <- arrange(res_sig, dir)


heatplot <- X[meta$SampleID,res_sig$feature_id] %>% t()


pheatmap(heatplot,
         scale = 'row',
         cluster_rows = F,
         cluster_cols = F
         )



#------------------------
# KEGG 富集分析结果
#------------------------

Pathwayres <- read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/AEE_PEE_metabolomics_OPLSDA_Wilcox/AEE KEGG enrichment/pathway_results.xlsx")

Pathwayres <- Pathwayres %>% subset(Raw.p < 0.1)

ggplot(Pathwayres, aes(-log10(Raw.p),
                       reorder(Pathway, Raw.p),
                       fill = Group))+
  geom_col(color = 'black', width = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        legend.position = 'inside',
        legend.position.inside = c(0.8, 0.8),
        axis.title.y = element_blank(),
        )+
  # facet_wrap(Group~., scale = 'free_y', ncol = 1)+
  scale_fill_manual(values = vol_color)



ggsave(file.path(outdir, "pathyway_enrichment_plot.pdf"), width = 5, height = 8)


#-------------------------------
# VFA与差异代谢物相关性分析
#--------------------------------

vfa_raw <- read.xlsx("./代谢组数据/VFA profiles.xlsx", sheet = 1)

# 期望的列名（按你表里的实际列名）
id_col <- "SampleID"
vfa_cols <- c("Acetic_acid","Propionic_acid","Isobutyric_acid",
              "Butyric_acid","Isovaleric_acid","Valeric_acid")


vfa_raw <- vfa_raw %>%
  rename_with(~trimws(.x)) %>%      # 清理列名空格（很重要）
  mutate(across(all_of(vfa_cols), as.numeric)) %>%
  group_by(.data[[id_col]]) %>%
  summarise(
    across(
      all_of(vfa_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


vfa <- left_join(ass0, vfa_raw, by = 'SampleID') %>% distinct()


met_mat <- heatplot <- X[meta$SampleID,res_sig$feature_id] %>% as.data.frame() %>% rownames_to_column(var = 'SampleID')

corr_mat <- inner_join(vfa,met_mat,  by = 'SampleID')


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(Hmisc)     # rcorr
  library(pheatmap)
})

# =========================
# 0) 你要相关的两类列
# =========================
# VFA / 发酵相关列（按你给的列名写死，避免误选）
vfa_cols <- c(
  # "TVFA","AP","ABP","BCVFA"
  "Acetic_acid","Propionic_acid","Isobutyric_acid",
  "Butyric_acid","Isovaleric_acid","Valeric_acid"
)

# 非靶代谢物列：pos_ / neg_ 开头（你也可以改成更严格的正则）
meta_cols <- grep("^(pos_|neg_)", colnames(corr_mat), value = TRUE)

# 检查列是否都存在
miss_vfa <- setdiff(vfa_cols, colnames(corr_mat))
if (length(miss_vfa) > 0) stop("corr_mat missing VFA columns: ", paste(miss_vfa, collapse = ", "))
if (length(meta_cols) == 0) stop("No non-target metabolite columns found with pattern ^(pos_|neg_)")

# 取子矩阵，并确保是数值
X_vfa  <- corr_mat[, vfa_cols, drop = FALSE] %>% mutate(across(everything(), as.numeric))
X_meta <- corr_mat[, meta_cols, drop = FALSE] %>% mutate(across(everything(), as.numeric))

# 合并后做 rcorr（它需要矩阵）
X_all <- cbind(X_vfa, X_meta)
X_all <- as.matrix(X_all)

# 可选：去掉全 NA 或方差为0的列（否则相关会报 NA）
is_bad <- apply(X_all, 2, function(x) all(is.na(x)) || sd(x, na.rm = TRUE) == 0)
if (any(is_bad)) {
  message("Dropping ", sum(is_bad), " columns that are all-NA or zero-variance.")
  X_all <- X_all[, !is_bad, drop = FALSE]
}

# 更新列集合（防止刚刚被删掉）
vfa_cols2  <- intersect(vfa_cols, colnames(X_all))
meta_cols2 <- intersect(meta_cols, colnames(X_all))

# =========================
# 1) 计算相关与P值（Spearman）
# =========================
rc <- Hmisc::rcorr(X_all, type = "spearman")
R  <- rc$r
P  <- rc$P

# 提取 VFA(行) × 非靶(列) 的子矩阵
R_sub <- R[vfa_cols2, meta_cols2, drop = FALSE]
P_sub <- P[vfa_cols2, meta_cols2, drop = FALSE]

# =========================
# 2) 多重校正（BH）得到 padj
# =========================
p_vec <- as.vector(P_sub)
# 注意：p_vec里会有NA（比如某些列有效样本太少），要保留NA位置
padj_vec <- rep(NA_real_, length(p_vec))
ok <- !is.na(p_vec)
padj_vec[ok] <- p.adjust(p_vec[ok], method = "BH")

padj_mat <- matrix(padj_vec, nrow = nrow(P_sub), ncol = ncol(P_sub),
                   dimnames = dimnames(P_sub))

# =========================
# 3) 星号矩阵：padj < 0.05 标 *
# =========================
stars <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat),
                dimnames = dimnames(padj_mat))
stars[!is.na(padj_mat) & padj_mat < 0.01] <- "*"

# （可选）如果你想分层星号：
# stars[!is.na(padj_mat) & padj_mat < 0.001] <- "***"
# stars[!is.na(padj_mat) & padj_mat >= 0.001 & padj_mat < 0.01] <- "**"
# stars[!is.na(padj_mat) & padj_mat >= 0.01 & padj_mat < 0.05] <- "*"

# =========================
# 4) 画热图：显示相关系数，并用 display_numbers 标星号
# =========================
# 可选：对列做聚类太多会很慢（你代谢物可能上千列），可以关闭列聚类或先筛选
# 例如只保留每个VFA相关性最高的前N个代谢物：
# N <- 80
# keep_meta <- unique(unlist(apply(abs(R_sub), 1, function(x) names(sort(x, decreasing = TRUE))[1:N])))
# R_sub <- R_sub[, keep_meta, drop = FALSE]
# stars <- stars[, keep_meta, drop = FALSE]

pheatmap::pheatmap(
  R_sub,
  display_numbers = stars,
  number_color = "black",
  fontsize_number = 9,
  cluster_rows = TRUE,
  cluster_cols = TRUE,     # 列太多可改 FALSE
  border_color = NA
)

# =========================
# 5) 输出结果表（可选）
# =========================
out_long <- as.data.frame(as.table(R_sub)) %>%
  rename(VFA = Var1, Metabolite = Var2, rho = Freq) %>%
  mutate(
    p = as.vector(P_sub),
    padj = as.vector(padj_mat),
    sig = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  ) %>%
  arrange(padj)

 write.csv(out_long, "VFA_vs_metabolites_spearman_BH.csv", row.names = FALSE)


# =========================
# 3.5) 强相关 + 显著性过滤
# =========================
rho_cut  <- 0.2
padj_cut <- 0.05

keep_mat <- (!is.na(R_sub)) &
  (abs(R_sub) >= rho_cut) &
  (!is.na(padj_mat)) &
  (padj_mat < padj_cut)

# 如果一个代谢物在所有 VFA 上都不满足条件，直接丢掉
keep_cols <- colSums(keep_mat) > 0
R_filt    <- R_sub[, keep_cols, drop = FALSE]
padj_filt <- padj_mat[, keep_cols, drop = FALSE]

# 如果一个 VFA 在所有代谢物上都没信号，也可以顺便丢掉（可选）
keep_rows <- rowSums(keep_mat[, keep_cols, drop = FALSE]) > 0
R_filt    <- R_filt[keep_rows, , drop = FALSE]
padj_filt <- padj_filt[keep_rows, , drop = FALSE]

# 安全检查
if (ncol(R_filt) == 0) stop("No metabolites pass |rho| >= 0.5 & padj < 0.05")


stars_filt <- matrix("", 
                     nrow = nrow(padj_filt), 
                     ncol = ncol(padj_filt),
                     dimnames = dimnames(padj_filt)
)

stars_filt[padj_filt < padj_cut] <- "*"

col_annotation <- res_sig


rownames(col_annotation) <- col_annotation$feature_id

pheatmap::pheatmap(
  R_filt,
  display_numbers = stars_filt,
  annotation_col = col_annotation[, c(32, 35), drop = FALSE],
  number_color = "black",
  fontsize_number = 9,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  show_colnames = F
)


out_strong <- as.data.frame(as.table(R_filt)) %>%
  rename(
    VFA = Var1,
    Metabolite = Var2,
    rho = Freq
  ) %>%
  mutate(
    padj = as.vector(padj_filt),
    direction = ifelse(rho > 0, "Positive", "Negative")
  ) %>%
  arrange(VFA, desc(abs(rho)))

write.csv(out_strong, "VFA_nonTarget_strong_correlations.csv", row.names = FALSE)


#-----------------------------------
# HMDB大类差异统计
#-----------------------------------

sig_sum <- table(res_sig$HMDB.Superclass, res_sig$Dir) %>% as.data.frame()


ggplot(sig_sum, aes(reorder(Var1, Freq), 
                    Freq, 
                    group = Var2, 
                    colour = Var2, fill = Var2))+
  geom_col(position = 'dodge', width = 0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    # axis.text.x = element_text(angle = 315, hjust = 0),
    panel.grid = element_blank(),
    axis.ticks.length = unit(2, "mm")
  ) +
  scale_fill_manual(values = c('AEE' = '#2e528f',
                                'PEE' = '#b41d23'))+
  scale_color_manual(values = c('AEE' = darken('#2e528f',0.5),
                                'PEE' = darken('#b41d23',0.5)))+
  coord_flip()


library(tidyverse)
library(scales)

# 1) 先决定排序（用总和排序，避免 AEE/PEE 各排各的）
ord <- sig_sum %>%
  group_by(Var1) %>%
  summarise(tot = sum(Freq, na.rm = TRUE), .groups = "drop") %>%
  arrange(tot) %>%
  pull(Var1)

sig_sum2 <- sig_sum %>%
  mutate(
    Var1 = factor(Var1, levels = ord),
    Var2 = factor(Var2, levels = c("AEE","PEE"))
  )

# 2) 画“圆环极坐标” + 并排柱
p <- ggplot(sig_sum2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(
    aes(color = Var2),
    position = position_dodge2(width = 0.75, preserve = "single"),
    width = 0.7,
    linewidth = 0.6
  ) +
  # 圆环：把 x 往右推，留中空
  coord_polar(start = -pi/2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(values = c(AEE = "#2e528f", PEE = "#b41d23")) +
  scale_color_manual(values = c(AEE = "#2e528f", PEE = "#b41d23")) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# 关键：中空圆环效果（控制“内径”）
p 



ggsave(file.path(outdir, "HMDBsuperclass of differential metabolites.pdf"), width = 8, height = 8)



sum_mat <- X[meta$SampleID,res_sig$feature_id] %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'ID')


sum_mat <- inner_join(met0[,c(1,2,23:25)], sum_mat, by = 'ID')








sum_mat_sum <- aggregate(sum_mat[,c(6:475)], by = list(HMDBsuperclass = sum_mat$HMDB.Superclass), sum)


sum_mat_sum <- sum_mat_sum %>% column_to_rownames(var = 'HMDBsuperclass') %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'SampleID')


sum_mat_sum <- right_join(meta,sum_mat_sum,by = 'SampleID')

sum_mat_sum_long <- pivot_longer(sum_mat_sum, 4:19, names_to = 'HMDBsuperclass', values_to = 'Abundance')


sum_mat_sum_long <- unite(sum_mat_sum_long,col = 'Group_HMDBsuperclass',HMDBsuperclass, Group,remove = F)

sum_mat_sum_long <- sum_mat_sum_long %>% subset(sum_mat_sum_long$HMDBsuperclass != 'Others')

ggplot(sum_mat_sum_long, aes(x = reorder(Group_HMDBsuperclass, desc(Abundance)), 
                             y = Abundance, 
                             color = Group, fill = Group)) +
  geom_point(
    size  = 1,            # 控制横杠长度
    alpha = 0.2
  ) +
  # 均值水平线段
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "errorbar",
    width = 0.3,       # 控制线段横向延伸长度
    linewidth = 1.5,        # 线段粗细
    show.legend = FALSE
  ) +
  # 标准差误差线
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8,
  )+
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 315, hjust = 0),
    panel.grid = element_blank(),
    axis.ticks.length = unit(2, "mm")
  ) +
  scale_color_manual(values = c('AEE' = '#2e528f',
                                'PEE' = '#b41d23'))

ggsave(file.path(outdir, "HMDBsuperclass of differential metabolites.pdf"), width = 8, height = 4)
  
  
  
  
library(ggplot2)

# 生成示例数据
set.seed(123)
data <- data.frame(value = rnorm(100))

# 绘制条形码图
ggplot(data, aes(x = value, y = 0)) +
  geom_rug(sides = "b", length = unit(0.1, "npc")) + # 底部条形码
  theme_minimal() +
  theme(axis.text.y = element_blank(), # 隐藏y轴
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(title = "Barcode Plot", x = "Value")

# 使用上面的数据，添加颜色分类
data$type <- ifelse(data$value > 0, "Positive", "Negative")

ggplot(data, aes(y = value, x = 0.5, fill = type)) +
  geom_tile(width = 0.01, height = 1) + # width控制条形宽度
  scale_fill_manual(values = c("Negative" = "red", "Positive" = "blue")) +
  theme_void() + # 极简主题
  theme(legend.position = "bottom") +
  labs(title = "Customized Barcode Plot")

  
  
# ---------- 7) WGCNA on differential metabolites ----------
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}
suppressPackageStartupMessages(library(WGCNA))

# WGCNA 推荐：关闭字符串因子 & 允许多线程
options(stringsAsFactors = FALSE)
allowWGCNAThreads()   # Windows 也能跑，只是线程可能受限

# 7.1 构建 WGCNA 输入矩阵：samples x diff-metabolites
diff_feats <- intersect(res_sig$feature_id, colnames(X))
if (length(diff_feats) < 30) {
  warning("差异代谢物数量 < 30，WGCNA 模块可能不稳定；建议放宽阈值或不加VIP过滤后再做。")
}

datExpr0 <- as.data.frame(X[, diff_feats, drop = FALSE])  # samples x features
# 可选：WGCNA 常用做法是对每个特征标准化（Z-score），尤其不同代谢物量纲差异大时
datExpr0 <- as.data.frame(scale(datExpr0))

# 清理坏样本/坏特征
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 7.2 准备 trait：AP（从 assignment 表里取）
# 你的 ass 里目前只有 SampleID/Rumenotype/Group；所以从 ass0 里尝试找 AP 列
ass0_ <- ass0 %>% mutate(SampleID = as.character(SampleID))

# 自动识别 AP 列名（你可以把候选列名按你的真实列名再加一些）
ap_candidates <- c("AP", "A_P", "Acetate_Propionate", "acetate_propionate",
                   "Acetate/Propionate", "Acetate.Propionate", "AP_ratio")

ap_col <- intersect(ap_candidates, colnames(ass0_))
if (length(ap_col) == 0) {
  stop("在 assignment 表里没找到 AP 列。请确认 ass0 里 AP 的列名，并把它加到 ap_candidates 里。")
}
ap_col <- ap_col[1]

trait_df <- ass0_ %>%
  transmute(
    SampleID = as.character(SampleID),
    AP = suppressWarnings(as.numeric(.data[[ap_col]]))
  ) %>%
  filter(!is.na(SampleID), !is.na(AP)) %>%
  distinct(SampleID, .keep_all = TRUE)

# 对齐样本顺序
common_samp <- intersect(rownames(datExpr0), trait_df$SampleID)
datExpr <- datExpr0[common_samp, , drop = FALSE]
datTraits <- trait_df %>%
  filter(SampleID %in% common_samp) %>%
  arrange(match(SampleID, common_samp)) %>%
  tibble::column_to_rownames("SampleID")

stopifnot(identical(rownames(datExpr), rownames(datTraits)))

# 7.3 选择 soft-threshold power
powers <- c(1:30)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "use='p'"
)

#拟合指数与 power 值散点图
ggplot(sft$fitIndices, aes(Power,SFT.R.sq, label=Power))+
  geom_line(linewidth=1,color='gray')+
  geom_point(size=2)+
  geom_text(color='red')+
  theme_bw()+
  geom_hline(yintercept = 0.8,color='red', linetype="dashed")+
  xlab ("Soft Threshold (power)")+
  ylab ("Scale Free Topology Model Fit")+
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank())

ggsave('./AEE_PEE_metabolomics_OPLSDA_Wilcox/WGCNA_SFT.pdf', width = 3, height = 3)


ggplot(sft$fitIndices, aes(Power,mean.k., label=Power))+
  geom_line(linewidth=1,color='gray')+
  geom_point(size=2)+
  geom_text(color='red')+
  theme_bw()+
  xlab ("Soft Threshold (power)")+
  ylab ("Mean Connectivity")+
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank())

ggsave('./AEE_PEE_metabolomics_OPLSDA_Wilcox/WGCNA_Mean Connectivity.pdf', width = 3, height = 3)


##构建网络
#上一步估计的最佳 power 值
powers_picked <- 20

#获得拓扑矩阵，详情 ?adjacency、?TOMsimilarity
adjacency <- adjacency(datExpr, power = powers_picked)
tom_sim <- TOMsimilarity(adjacency)
rownames(tom_sim) <- rownames(adjacency)
colnames(tom_sim) <- colnames(adjacency)
tom_sim[1:6,1:6]

#输出拓扑矩阵
#write.table(tom_sim, 'TOMsimilarity.txt', sep = '\t', col.names = NA, quote = FALSE)

#查看此时的网络的无标度拓扑特征
k <- softConnectivity(datE = datExpr, power = powers_picked)
par(mfrow = c(1, 2))
hist(k)  #会呈现一种幂律分布的状态
scaleFreePlot(k, main = 'Check Scale free topology')

##共表达模块划分
#相异度 = 1 - 相似度
tom_dis  <- 1 - tom_sim

#层次聚类树，使用中值的非权重成对组法的平均聚合聚类
geneTree <- hclust(as.dist(tom_dis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)

#使用动态剪切树挖掘模块，详情 ?cutreeDynamic
minModuleSize <- 10  #模块基因数目
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)

#模块颜色指代
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
genecolors <- ifelse(res_sig$Dir=='AEE', '#2e528f','#b41d23')

plotDendroAndColors(geneTree, cbind(dynamicColors,genecolors), c('Dynamic Tree Cut','Group'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
                    main = 'Gene dendrogram and module colors')

#基因表达聚类树和共表达拓扑热图，详情 ?TOMplot
plot_sim <- -(1-tom_sim)
#plot_sim <- log(tom_sim)
diag(plot_sim) <- NA
TOMplot(plot_sim, geneTree,  
        Colors = dynamicColors, 
        ColorsLeft = genecolors,
        main = 'Network heatmap') 

##模块特征基因
#计算基因表达矩阵中模块的特征基因（第一主成分），详情 ?moduleEigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]

#输出模块特征基因矩阵
#write.table(MEs, 'moduleEigengenes.txt', sep = '\t', col.names = NA, quote = FALSE)

##共表达模块的进一步聚类
#通过模块特征基因计算模块间相关性，表征模块间相似度
ME_cor <- cor(MEs)
ME_cor[1:6,1:6]

#绘制聚类树观察
METree <- hclust(as.dist(1-ME_cor), method = 'average')
plot(METree, main = 'Clustering of module eigengenes')

#探索性分析，观察模块间的相似性
#height 值可代表模块间的相异度，并确定一个合适的阈值作为剪切高度
#以便为低相异度（高相似度）的模块合并提供依据
abline(h = 0.1, col = 'blue')
abline(h = 0.25, col = 'red')

#模块特征基因聚类树热图
plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90, 
                      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2))

#相似模块合并，以 0.25 作为合并阈值（剪切高度），在此高度下的模块将合并
#近似理解为相关程度高于 0.75 的模块将合并到一起
merge_module <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.1, verbose = 3)
mergedColors <- merge_module$colors

table(mergedColors)

#基因表达和模块聚类树
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors,genecolors), c('Dynamic Tree Cut', 'Merged dynamic','Group'),
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)


#使用上一步新组合的共表达模块的结果
module <- merge_module$newMEs

ass_AP <- ass0 %>% column_to_rownames(var = 'SampleID')
ass_AP <- ass_AP[rownames(module),]

module <- cbind(ass_AP[,4:7],module)

# 计算各个模块之间的相关性
# 计算相关性  
correlation_ms <- cor(module, method = 'spearman')
pvalue_ms <- corrplot::cor.mtest(module, method = 'spearman')$p


# 假设 pvalue_ms 是你的原始 p 值矩阵
p_vector <- pvalue_ms[lower.tri(pvalue_ms)]  # 提取下三角（不含对角线）
adjusted_p <- p.adjust(p_vector, method = "fdr")  # BH 方法矫正

# 构建矫正后的矩阵
adjusted_matrix <- pvalue_ms
adjusted_matrix[lower.tri(adjusted_matrix)] <- adjusted_p
adjusted_matrix[upper.tri(adjusted_matrix)] <- t(adjusted_matrix)[upper.tri(adjusted_matrix)]

# 保留行列名
rownames(adjusted_matrix) <- rownames(pvalue_ms)
colnames(adjusted_matrix) <- colnames(pvalue_ms)

# 查看结果
head(adjusted_matrix)

# 假设 adjusted_matrix 是矫正后的数值矩阵
# 创建新矩阵存储星号标记，保持原始数值矩阵不变
symbol_matrix <- matrix("", nrow = 23, ncol = 23)
rownames(symbol_matrix) <- rownames(adjusted_matrix)
colnames(symbol_matrix) <- colnames(adjusted_matrix)

# 定义映射规则
for (i in 1:23) {
  for (j in 1:23) {
    p <- adjusted_matrix[i, j]
    if (i == j) {
      symbol_matrix[i, j] <- ""  # 对角线留空
    } else if (p <= 0.001) {
      symbol_matrix[i, j] <- "***"
    } else if (p <= 0.01) {
      symbol_matrix[i, j] <- "**"
    } else if (p <= 0.05) {
      symbol_matrix[i, j] <- "*"
    } else {
      symbol_matrix[i, j] <- "NS"
    }
  }
}

diag(correlation_ms) <- NA

library(pheatmap)

# 1. 创建注释数据框（假设行和列注释为模块名称）
annot_df <- data.frame(
  Module = rownames(correlation_ms)
)
rownames(annot_df) <- rownames(correlation_ms)

annot_color <- setNames(
  c(
    rep("grey", 4),
    sub("^ME", "", c("MEdarkgrey","MEroyalblue","MEdarkturquoise","MEgrey60","MEpink",
                     "MEblack","MEsalmon","MEmagenta","MEmidnightblue","MElightyellow",
                     "MElightcyan","MEbrown","MEred","MEcyan","MElightgreen",
                     "MEdarkred","MEdarkgreen","MEtan","MEgrey"))
  ),
  c("TVFA","AP","ABP","BCVFA",
    "MEdarkgrey","MEroyalblue","MEdarkturquoise","MEgrey60","MEpink",
    "MEblack","MEsalmon","MEmagenta","MEmidnightblue","MElightyellow",
    "MElightcyan","MEbrown","MEred","MEcyan","MElightgreen",
    "MEdarkred","MEdarkgreen","MEtan","MEgrey")
)


# 3. 绘制热图
pheatmap(
  mat = correlation_ms[-12, -12],
  display_numbers = symbol_matrix[-12, -12],
  annotation_row = annot_df,    # 行注释数据框
  annotation_col = annot_df,    # 列注释数据框
  annotation_colors = list(
    Module = annot_color        # 颜色映射到分类变量 "Module"
  )
)

# 经验选择：优先使 scale-free fit >= 0.8 的最小 power；否则取使其“拐点”的 power
fit_df <- sft$fitIndices
if (is.null(fit_df)) stop("pickSoftThreshold failed: sft$fitIndices is NULL.")
candidate <- fit_df$Power[fit_df$SFT.R.sq >= 0.80]
softPower <- if (length(candidate) > 0) min(candidate) else fit_df$Power[which.max(fit_df$SFT.R.sq)]
message("Chosen softPower = ", softPower)

# 7.4 构建网络 + 模块识别（blockwiseModules 更稳健）
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 20,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs0 <- net$MEs
MEs <- orderMEs(MEs0)

# 7.5 模块-性状相关：ME vs AP
# 用 bicor 更稳健；p 值用 corPvalueStudent（近似）
ME_AP_cor <- bicor(MEs, datTraits$AP, use = "p")
ME_AP_p   <- corPvalueStudent(ME_AP_cor, nSamples = nrow(datExpr))

ME_AP_tab <- data.frame(
  Module = colnames(MEs),
  cor_AP = as.numeric(ME_AP_cor),
  p_AP   = as.numeric(ME_AP_p),
  padj_AP = p.adjust(as.numeric(ME_AP_p), method = "BH"),
  stringsAsFactors = FALSE
) %>% arrange(p_AP)

# 显著模块（你可以改阈值）
sig_mods <- ME_AP_tab %>% filter(padj_AP < 0.05)
print(sig_mods)

# 7.6 输出：模块-AP 相关热图
pdf(file.path(outdir, "WGCNA_ModuleTrait_ME_vs_AP.pdf"), width = 6, height = 7)
textMatrix <- paste0(signif(ME_AP_cor, 2), "\n(",
                     signif(ME_AP_p, 2), ")")
dim(textMatrix) <- dim(ME_AP_cor)
labeledHeatmap(
  Matrix = ME_AP_cor,
  xLabels = "AP",
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module eigengenes vs AP"
)
dev.off()

# 7.7 代谢物归属（feature -> module）+ 模块成员度（kME）+ hub metabolites
kME <- as.data.frame(bicor(datExpr, MEs, use = "p"))
colnames(kME) <- paste0("kME_", colnames(MEs))

feat_anno <- res_sig %>%
  select(feature_id, everything()) %>%
  distinct(feature_id, .keep_all = TRUE)

module_assign <- data.frame(
  feature_id = colnames(datExpr),
  moduleColor = moduleColors,
  moduleLabel = moduleLabels,
  stringsAsFactors = FALSE
) %>%
  left_join(feat_anno, by = "feature_id") %>%
  bind_cols(kME[match(colnames(datExpr), rownames(kME)), , drop = FALSE])

# 每个模块的 hub：按对应模块的 kME 降序取前 20
get_top_hubs <- function(df, module_color, top_n = 20) {
  mod_name <- paste0("ME", module_color)
  kme_col  <- paste0("kME_", mod_name)
  if (!kme_col %in% colnames(df)) return(NULL)
  df %>%
    filter(moduleColor == module_color) %>%
    arrange(desc(.data[[kme_col]])) %>%
    slice_head(n = top_n) %>%
    select(feature_id, moduleColor, all_of(kme_col),
           any_of(c("VIP", "p_wilcox", "padj", "Metabolite", "mz", "rt", "Mode", "Dir")))
}

# 导出全量归属表
write.csv(module_assign, file.path(outdir, "WGCNA_feature_module_membership.csv"), row.names = FALSE)

# 导出显著模块列表 + 每个显著模块 hub metabolites
write.csv(ME_AP_tab, file.path(outdir, "WGCNA_ME_AP_cor_table.csv"), row.names = FALSE)

if (nrow(sig_mods) > 0) {
  for (m in sig_mods$Module) {
    # m 形如 "MEturquoise"；取 color 名称
    mod_color <- sub("^ME", "", m)
    hubs <- get_top_hubs(module_assign, mod_color, top_n = 30)
    if (!is.null(hubs)) {
      write.csv(hubs, file.path(outdir, paste0("WGCNA_Hubs_", m, "_top30.csv")), row.names = FALSE)
    }
  }
}

#-----------------------------
#差异代谢物人工标注，以及富集分析
#-----------------------------

res_sig_label <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 2-AEE PEE differential metabolites.xlsx")
res_all <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/AEE_PEE_metabolomics_OPLSDA_Wilcox/results_all_metabolites.tsv")



# (3) Volcano: delta_med vs -log10(p)
res_plot <- res_all %>%
  mutate(
    neglog10p = -log10(p_wilcox),
    Sig = ifelse(p_wilcox < p_cut & VIP >1 , "p<0.05", "NS"),
    Dir = ifelse(p_wilcox < p_cut & VIP >1 & median_AEE > median_PEE, 'AEE',
                 ifelse(p_wilcox < p_cut & VIP >1 & median_AEE < median_PEE, 'PEE','N.S.'
                 ))
  )

res_plot <- left_join(res_plot, res_sig_label[,c(1,37)], by='feature_id')

ferment_cols <- c(
  "Peptide- and amino acid–related fermentation products" = "#2F5597",  # deep blue
  "Carbohydrate fermentation intermediates"               = "#1B998B",  # teal-green
  "Lipid-related / fatty acid–associated metabolites"     = "#C97C5D",  # muted amber/terracotta
  # "Reducing equivalents sinks / redox-related"            = "#6A4C93",   # elegant purple
  'Sig'='#969696', 
  'N.S.'= '#efefef'
)

res_plot$fermentation_class <- ifelse(!is.na(res_plot$fermentation_class), res_plot$fermentation_class,
                                      ifelse(is.na(res_plot$fermentation_class) & 
                                               res_plot$Dir != 'N.S.', 'Sig', 'N.S.'
                                             )
                                      )


p_vol <- ggplot(res_plot, aes(x = log2fc, y = -log10(p_wilcox), colour = fermentation_class)) +
  geom_point( size = 2) +
  theme_bw(base_size = 12) +
  # labs( x = "delta median (AEE - PEE)", y = "-log10(Wilcoxon p)") +
  geom_hline(yintercept = -log10(p_cut), linetype = 2) +
  scale_color_manual(values = ferment_cols)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        axis.ticks.length = unit(2,'mm'))

p_vol

table(res_plot$fermentation_class)

ggsave(file.path(outdir, "Volcano_fermentation_class.pdf"), p_vol, width = 5, height = 5)



#------------------
#基于3类中间代谢物绘制三元相图
#------------------

# 合并到 res_plot（res_plot 的 feature_id 类似 pos_1 / pos_976）
res_gsea <- res_plot %>%
  left_join(fer_class_df, by = c("feature_id" = "ID"))

# 只保留三类（你要的三类）
keep_classes <- c(
  "Carbohydrate fermentation intermediates",
  "Peptide- and amino acid–related fermentation products",
  "Lipid-related / fatty acid–associated metabolites"
)

res_gsea <- res_plot %>% filter(fermentation_class %in% keep_classes)

pie_plot <- table(res_gsea$fermentation_class, res_gsea$Dir) %>% as.data.frame()

group_cols <- c(
  AEE = "#2e528f",  # 你之前常用的红
  PEE = "#b41d23"   # 青蓝
)


library(tidyverse)

pie_df <- pie_plot %>%
  group_by(Var1) %>%
  mutate(
    prop = Freq / sum(Freq),
    label = sprintf("%.1f%%", prop * 100)
  ) %>%
  ungroup()


plot_one_pie <- function(df, title_text) {
  ggplot(df, aes(x = "", y = Freq, fill = Var2)) +
    geom_col(width = 1, color = "white", size = 0.6) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = group_cols) +
    theme_void(base_size = 12) +
    labs(title = title_text) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.title = element_blank(),
      legend.position = "right"
    )
}

p_carb <- pie_df %>%
  filter(Var1 == "Carbohydrate fermentation intermediates") %>%
  plot_one_pie("Carbohydrate-oriented metabolites")

p_carb

ggsave(file.path(outdir, "p_carb.pdf"), p_carb, width = 5, height = 5)

p_lipid <- pie_df %>%
  filter(Var1 == "Lipid-related / fatty acid–associated metabolites") %>%
  plot_one_pie("Lipid-oriented metabolites")

p_lipid

ggsave(file.path(outdir, "p_lipid.pdf"), p_lipid, width = 5, height = 5)

p_aa <- pie_df %>%
  filter(Var1 == "Peptide- and amino acid–related fermentation products") %>%
  plot_one_pie("Amino acid–oriented metabolites")

p_aa

ggsave(file.path(outdir, "p_aa.pdf"), p_aa, width = 5, height = 5)
