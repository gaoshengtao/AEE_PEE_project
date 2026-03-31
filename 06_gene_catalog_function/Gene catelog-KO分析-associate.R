# =========================
# KO-based functional analysis (FINAL):
#   TPM -> CLR
#   limma: ~ scale(AP)            -> beta_AP_CLR
#   limma: ~ Rumenotype           -> beta_group_CLR (Acetate - Propionate)
#   Consistency: correlation + quadrant plot (same-sign fraction)
#   limma: ~ Rumenotype * scale(AP) -> interaction slope difference
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(limma)
  library(ggplot2)
  library(scales)
  library(readr)
  library(tibble)
})

# --------- 0) Inputs (EDIT) ----------
ko_tsv   <- "./微生物/genecatelog/06_aggregate_tpm_by_eggnog/TPM_by_KEGG_ko.tsv"   # KO x Sample matrix (tsv)
meta_xlsx <- "./FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx" # meta: SampleID, Rumenotype, AP
outdir   <- "./微生物/KO_AEE_PEE_CLR_LIMMA"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# metadata columns
sample_col <- "SampleID"
group_col  <- "Rumenotype"   # levels: Acetate_type / Propionate_type
ap_col     <- "AP"

# filtering & CLR params
min_prevalence <- 0.10  # keep KO present in >=10% samples (TPM > 0)
pseudo_count   <- 0.5   # added before log for CLR
set.seed(1)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clr_transform <- function(mat, pseudo = 0.5) {
  # mat: features x samples, non-negative
  mat2 <- mat + pseudo
  # log
  logm <- log(mat2)
  # CLR: subtract per-sample mean(log)
  sweep(logm, 2, colMeans(logm, na.rm = TRUE), FUN = "-")
}

safe_write <- function(df, path) {
  readr::write_tsv(df, path)
  message("[write] ", path)
}

# ------------------------------------------------------------
# 1) Read metadata
# ------------------------------------------------------------
meta <- openxlsx::read.xlsx(meta_xlsx) %>%
  mutate(
    !!sample_col := as.character(.data[[sample_col]]),
    !!group_col  := as.factor(.data[[group_col]]),
    !!ap_col     := as.numeric(.data[[ap_col]])
  )

stopifnot(sample_col %in% names(meta), group_col %in% names(meta), ap_col %in% names(meta))

# enforce group level order (important for sign)
meta[[group_col]] <- factor(meta[[group_col]], levels = c("Acetate_type", "Propionate_type"))
message("[meta] group levels: ", paste(levels(meta[[group_col]]), collapse = ", "))

# ------------------------------------------------------------
# 2) Read KO TPM matrix
# ------------------------------------------------------------
ko <- fread(ko_tsv)
names(ko)[1] <- "KO"

sample_cols <- intersect(names(ko), meta[[sample_col]])
stopifnot(length(sample_cols) >= 4)

ko_mat <- ko %>%
  select(KO, all_of(sample_cols)) %>%
  as.data.frame()

rownames(ko_mat) <- ko_mat$KO
ko_mat$KO <- NULL

# align meta to matrix columns
meta2 <- meta %>%
  filter(.data[[sample_col]] %in% colnames(ko_mat)) %>%
  arrange(match(.data[[sample_col]], colnames(ko_mat)))

ko_mat <- ko_mat[, meta2[[sample_col]], drop = FALSE]

stopifnot(all(colnames(ko_mat) == meta2[[sample_col]]))

# ------------------------------------------------------------
# 3) Filter low-prevalence KOs (based on TPM > 0)
# ------------------------------------------------------------
prev <- rowMeans(ko_mat > 0, na.rm = TRUE)
keep <- prev >= min_prevalence
ko_f <- ko_mat[keep, , drop = FALSE]

message(sprintf("[filter] KOs: %d (raw) -> %d (prevalence >= %.0f%%)",
                nrow(ko_mat), nrow(ko_f), 100 * min_prevalence))

# ------------------------------------------------------------
# 4) TPM -> CLR
# ------------------------------------------------------------
clr_mat <- clr_transform(as.matrix(ko_f), pseudo = pseudo_count)

# Remove rows with zero variance (can cause NA coef)
ko_var <- apply(clr_mat, 1, var, na.rm = TRUE)
clr_mat <- clr_mat[ko_var > 0, , drop = FALSE]

# sanity
stopifnot(all(is.finite(clr_mat)))
message("[CLR] matrix dim: ", nrow(clr_mat), " x ", ncol(clr_mat))

# Prepare meta rownames for limma
meta2 <- meta2 %>% mutate(SampleID = as.character(.data[[sample_col]]))
rownames(meta2) <- meta2$SampleID
stopifnot(all(colnames(clr_mat) == rownames(meta2)))

# ------------------------------------------------------------
# 5) limma model 1: CLR ~ scale(AP)   -> beta_AP_CLR
# ------------------------------------------------------------
design_ap <- model.matrix(~ scale(AP), data = meta2)
fit_ap <- lmFit(clr_mat, design_ap)
fit_ap <- eBayes(fit_ap)

coef_ap <- "scale(AP)"
res_ap <- topTable(fit_ap, coef = coef_ap, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("KO") %>%
  as_tibble() %>%
  transmute(
    KO,
    beta_AP_CLR = logFC,
    pvalue_AP   = P.Value,
    padj_AP     = adj.P.Val
  ) %>%
  arrange(padj_AP, pvalue_AP)

res_ap$color <- ifelse(res_ap$beta_AP_CLR > 0 & res_ap$pvalue_AP < 0.05, 'AP_pos',
                       ifelse(res_ap$beta_AP_CLR < 0 & res_ap$pvalue_AP < 0.05, 'AP_neg', 'N.S.'))

safe_write(res_ap, file.path(outdir, "KO_betaAP_CLR_limma.tsv"))

vol_color <- c('AP_pos'='#2e528f', 'AP_neg'='#b41d23', 'N.S.'= '#efefef')

# volcano-like scatter
p_ap <- ggplot(res_ap, aes(x = beta_AP_CLR, y = -log10(pvalue_AP), color = color)) +
  geom_point(size = 1) +
  # ggrepel::geom_text_repel(
  #   data = label_df,
  #   aes(label = mo),
  #   size = 3,
  #   max.overlaps = Inf
  # ) +
  labs(
    x = "Coefficient",
    y = '-Log10(p)'
  ) +
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm')
  )+
  scale_color_manual(values = vol_color)

p_ap

ggsave(file.path(outdir, "KO_betaAP_CLR_vol_scatter.pdf"), p_ap, width = 4, height = 4)

table(res_ap$color)

# ------------------------------------------------------------
# 6) limma model 2: CLR ~ Rumenotype  -> beta_group_CLR (Acetate - Propionate)
# ------------------------------------------------------------
meta2[[group_col]] <- factor(meta2[[group_col]], levels = c("Acetate_type", "Propionate_type"))

design_g <- model.matrix(~ Rumenotype, data = meta2)
fit_g <- lmFit(clr_mat, design_g)
fit_g <- eBayes(fit_g)

# coefficient represents (Propionate - Acetate); we want (Acetate - Propionate) => negate
coef_g <- "RumenotypePropionate_type"
stopifnot(coef_g %in% colnames(design_g))

tmp_g <- topTable(fit_g, coef = coef_g, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("KO") %>%
  as_tibble()

res_g <- tmp_g %>%
  transmute(
    KO,
    beta_group_CLR = -logFC,     # Acetate - Propionate
    pvalue_group   = P.Value,
    padj_group     = adj.P.Val
  ) %>%
  arrange(padj_group, pvalue_group)

res_g$color <- ifelse(res_g$beta_group_CLR > 0 & res_g$pvalue_group < 0.05, 'AP_pos',
                       ifelse(res_g$beta_group_CLR < 0 & res_g$pvalue_group < 0.05, 'AP_neg', 'N.S.'))

safe_write(res_g, file.path(outdir, "KO_betaGroup_CLR_AminusP_limma.tsv"))

table(res_g$color)

# volcano-like scatter
p_g <- ggplot(res_g, aes(x = beta_group_CLR, y = -log10(pvalue_group), color = color)) +
  geom_point(size = 1) +
  # ggrepel::geom_text_repel(
  #   data = label_df,
  #   aes(label = mo),
  #   size = 3,
  #   max.overlaps = Inf
  # ) +
  labs(
    x = "Coefficient",
    y = '-Log10(padj)'
  ) +
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm')
  )+
  scale_color_manual(values = vol_color)

p_g

ggsave(file.path(outdir, "KO_betagroup_CLR_vol_scatter.pdf"), p_g, width = 4, height = 4)

table(res_g$color)

# ------------------------------------------------------------
# 7) Consistency: correlation + quadrant plot (same-sign fraction)
# ------------------------------------------------------------
df_merge <- res_ap %>%
  inner_join(res_g, by = "KO") %>%
  mutate(
    consistent = (beta_AP_CLR * beta_group_CLR) > 0,
    cons_lab = ifelse(consistent, "Consistent", "Inconsistent")
  )

df_merge$cons_lab <- ifelse(df_merge$color.x == 'N.S.' & df_merge$color.y == 'N.S.', 'N.S.', df_merge$cons_lab)

safe_write(df_merge, file.path(outdir, "KO_betaAP_vs_betaGroup_CLR_merged.tsv"))

cons_summary <- df_merge %>%
  summarise(
    n = n(),
    prop_consistent = mean(consistent, na.rm = TRUE),
    spearman = cor(beta_AP_CLR, beta_group_CLR, method = "spearman", use = "pairwise.complete.obs"),
    pearson  = cor(beta_AP_CLR, beta_group_CLR, method = "pearson",  use = "pairwise.complete.obs")
  )
safe_write(cons_summary, file.path(outdir, "consistency_summary.tsv"))

COL_CONS <- c("Consistent" = "#D62728", "Inconsistent" = "#1F77B4", "N.S."='#efefef')

p_quad <- ggplot(df_merge, aes(x = beta_group_CLR, y = beta_AP_CLR, color = cons_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_point(alpha = 0.65, size = 1.05) +
  scale_color_manual(values = COL_CONS) +
  theme_bw(base_size = 11) +
  labs(
    x = "beta_group_CLR (Acetate − Propionate)",
    y = "beta_AP_CLR (effect of AP)",
    title = sprintf("KO consistency (CLR+limma): same-sign=%.1f%% | Spearman=%.3f",
                    100 * cons_summary$prop_consistent, cons_summary$spearman)
  )+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm')
  )

p_quad

ggsave(file.path(outdir, "KO_betaAP_vs_betaGroup_CLR_quadrant.pdf"),
       p_quad, width = 7, height = 5.0)

table(df_merge$cons_lab)


# 

module_name <- read.table("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/module name.txt", header = T, sep = '\t')

kegg_module_hierarchy <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/kegg_module_hierarchy_L1L2L3.tsv")

infer_module <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/06_infer_kegg_names_final/modules_bestname.tsv")

infer_module <- left_join(infer_module, kegg_module_hierarchy, by = c('Best_KEGG_Module' = 'module_id'))

ko_module <- read.table("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/module_ko_map.txt", header = T)

ko_pathway <- read.table("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/pathway_ko_map.txt", header = T)

infer_pathway <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/06_infer_kegg_names_final/pathways_bestname.tsv")

pathway_hierarchy <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/kegg_pathway_hierarchy_L1L2L3.tsv")

pathway_hierarchy$pathway_id <- gsub('^map','ko',pathway_hierarchy$pathway_id)

ko_pathway <- ko_pathway %>% subset(!Pathway %>% str_detect('map'))


consitent_ko <- df_merge %>% subset(cons_lab == 'Consistent')

consitent_ko$dir <- ifelse(consitent_ko$beta_AP_CLR > 0, 'AEE','PEE')

consitent_ko$KO <- gsub('^ko:','',consitent_ko$KO)


consitent_ko <- left_join(consitent_ko, ko_module, by='KO')



consitent_ko <- left_join(consitent_ko, ko_pathway, by='KO')

consitent_ko_pathway <- consitent_ko %>%
  select(1,12,14) %>%
  mutate(across(where(is.character), stringr::str_trim)) %>%
  distinct()


consitent_ko_pathway <- left_join(consitent_ko_pathway, infer_pathway, by = c('Pathway' = 'Best_KEGG_Pathway'))

consitent_ko_pathway <- left_join(consitent_ko_pathway, pathway_hierarchy, by = c('Pathway' = 'pathway_id'))

consitent_ko_pathway <- consitent_ko_pathway[,c(1:3,16:18)]

consitent_ko_pathway <- consitent_ko_pathway %>% distinct()
  

consitent_ko_pathway_level2 <- consitent_ko_pathway %>%
  select(1,2,4,5) %>%
  distinct()

consitent_ko_pathway_level2$level2[is.na(consitent_ko_pathway_level2$level2)] <- 'Others'

consitent_ko_pathway_level2_sum <- table(consitent_ko_pathway_level2$KO, consitent_ko_pathway_level2$level2) %>% as.data.frame() %>% 
  pivot_wider(names_from = 'Var2', values_from = 'Freq') %>% column_to_rownames(var = 'Var1')

top_level2 <- table(consitent_ko_pathway$level2) %>% sort(decreasing = T)

top_level2 <- top_level2[1:30] %>% names()

consitent_ko_pathway_level2_sum2 <- consitent_ko_pathway_level2_sum %>%
  mutate(Others = rowSums(across(-all_of(top_level2)), na.rm = TRUE)) %>%
  select(all_of(top_level2), Others)


col_annotation <- consitent_ko_pathway[,1:2] %>% distinct() %>% arrange(dir) %>% as.data.frame()

rownames(col_annotation) <- col_annotation$KO


rownames(ko_mat) <- gsub('^ko:','',rownames(ko_mat))


ko_mat_consitent <- ko_mat[col_annotation$KO,]



library(pheatmap)

# 1) 数值矩阵
mat <- as.matrix(ko_mat_consitent)
mode(mat) <- "numeric"

# 2) 你说：col_annotation 注释行（对应 mat 的 rownames）
#    meta2 注释列（对应 mat 的 colnames）

# ---- 行注释：取 col_annotation 第2列，但一定要保持 data.frame + 行名对齐 ----
ann_row <- col_annotation[rownames(mat), 2, drop = FALSE]
colnames(ann_row) <- "Row_anno"  # 可选：给注释列起个名字，legend更清晰

# ---- 列注释：取 meta2 第4:8列，保持 data.frame + 行名对齐 ----
ann_col <- meta2[colnames(mat), 4:8, drop = FALSE]

ann_col <- arrange(ann_col, Rumenotype)

mat <- mat[,rownames(ann_col)]

# 3) 画图（行标准化）
pheatmap(
  mat,
  scale = "row",
  annotation_row = ann_row,
  annotation_col = ann_col,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
)




# 1) 行标准化
mat_z <- t(scale(t(mat)))
mat_z[!is.finite(mat_z)] <- 0

# 2) 截断极值（避免被极端值控制）
cap <- 2
mat_z2 <- pmax(pmin(mat_z, cap), -cap)

# 3) 自定义 breaks & colors（保证色阶均匀）
bk <- seq(-cap, cap, length.out = 101)
cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# 4) 注释（按你方向：行注释 col_annotation；列注释 meta2）
ann_row <- col_annotation[rownames(mat_z2), 2, drop = FALSE]
ann_col <- meta2[colnames(mat_z2), 4:8, drop = FALSE]

pheatmap(
  mat_z2,
  color = cols,
  breaks = bk,
  annotation_row = ann_row,
  annotation_col = ann_col,
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  border_color = NA
)



























# --------- 3) PCoA (CLR + Euclidean = Aitchison distance) ----------
# CLR transform on samples (compositional)
# NOTE: if ko_f are counts/abundances, this is fine for ordination
x <- t(as.matrix(ko_f))  # samples x KOs
x <- x + pseudo_count

# closure (to relative abund) then CLR
x_rel <- sweep(x, 1, rowSums(x), "/")
clr <- function(v) log(v) - mean(log(v))
x_clr <- t(apply(x_rel, 1, clr))  # samples x KOs

dist_aitch <- dist(x_clr, method = "euclidean")
pcoa <- cmdscale(dist_aitch, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  SampleID = rownames(pcoa$points),
  PCoA1 = pcoa$points[,1],
  PCoA2 = pcoa$points[,2]
) %>%
  left_join(meta2, by = c("SampleID" = sample_col))

# adonis2 + betadisper (建议一起报，避免“离散度差异”被质疑)
adon <- adonis2(dist_aitch ~ Rumenotype, data = pcoa_df, permutations = 999)
adon
bd   <- betadisper(dist_aitch, group = pcoa_df[[group_col]])
bd_p <- permutest(bd, permutations = 999)

fwrite(as.data.frame(adon), file.path(outdir, "adonis2_Aitchison.tsv"), sep = "\t")
fwrite(data.frame(betadisper_F = bd_p$tab[1,"F"], betadisper_p = bd_p$tab[1,"Pr(>F)"]),
       file.path(outdir, "betadisper.tsv"), sep = "\t")

# PCoA plot
var1 <- round(100 * pcoa$eig[1] / sum(pcoa$eig), 1)
var2 <- round(100 * pcoa$eig[2] / sum(pcoa$eig), 1)

p_pcoa <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = .data[[group_col]])
) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(
    aes(fill = .data[[group_col]]),
    type = "t",
    level = 0.95,
    geom = "polygon",
    alpha = 0
  ) +
  labs(
    # subtitle = sprintf("adonis2: R2=%.3f, p=%.3g | betadisper p=%.3g",adon$R2[1], adon$`Pr(>F)`[1], bd_p$tab[1, "Pr(>F)"]),
    x = paste0("PCoA1 (", var1, "%)"),
    y = paste0("PCoA2 (", var2, "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.8,0.8),
    # axis.text = element_blank(),
    axis.ticks.length = unit(2,'mm'),
    panel.grid = element_blank()
  )+
  scale_color_manual(values = c('Acetate_type'='#2e528f',
                                'Propionate_type'='#b41d23'))

p_pcoa


ggsave(file.path(outdir, "PCoA_Aitchison_KO.pdf"), p_pcoa, width = 4, height = 4)


#---------------------------------------------
# 差异KEGG ko的富集分析
#-------------------------------

# =========================
# (NEW) Hypergeometric enrichment on KEGG pathways
# =========================
P_CUT <- 0.05
ENRICH_FDR <- 0.1
MIN_K <- 3

df_merge$KO <- sub('^ko:','',df_merge$KO)

# universe
universe_ko <- df_merge$KO

# nominal diff KO
sig_ko <- df_merge %>%
  dplyr::mutate(
    KO  = sub("^ko:", "", KO),
    dir = dplyr::case_when(
      beta_AP_CLR*beta_group_CLR > 0 & beta_AP_CLR > 0 & (pvalue_AP < P_CUT| pvalue_group< P_CUT)  ~ "AEE",
      beta_AP_CLR*beta_group_CLR > 0 & beta_AP_CLR < 0 & (pvalue_AP < P_CUT| pvalue_group< P_CUT)  ~ "PEE",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(dir), pvalue_AP < P_CUT | pvalue_group < P_CUT) %>%
  dplyr::select(KO, dir) %>%
  dplyr::distinct()

# KO -> Pathway (map only) within universe
ko2path_all <- ko_pathway %>%
  dplyr::mutate(
    KO = sub("^ko:", "", KO),
    Pathway = as.character(Pathway)
  ) %>%
  dplyr::filter(stringr::str_detect(Pathway, "^map\\d{5}$")) %>%
  dplyr::filter(KO %in% universe_ko) %>%
  dplyr::distinct(KO, Pathway)

universe_mappable <- unique(ko2path_all$KO)
N <- length(universe_mappable)

path_bg <- ko2path_all %>%
  dplyr::count(Pathway, name = "M")

run_hyper <- function(dir_label) {
  query_ko <- sig_ko %>%
    dplyr::filter(dir == dir_label) %>%
    dplyr::pull(KO) %>%
    unique()
  
  query_ko <- intersect(query_ko, universe_mappable)
  n <- length(query_ko)
  
  hit <- ko2path_all %>%
    dplyr::filter(KO %in% query_ko) %>%
    dplyr::count(Pathway, name = "k")
  
  path_bg %>%
    dplyr::left_join(hit, by = "Pathway") %>%
    dplyr::mutate(
      k = tidyr::replace_na(k, 0L),
      n = n,
      N = N,
      p = stats::phyper(q = k - 1, m = M, n = N - M, k = n, lower.tail = FALSE),
      dir = dir_label
    ) %>%
    dplyr::mutate(padj = p.adjust(p, method = "BH")) %>%
    dplyr::arrange(p)
}

enrich_all <- dplyr::bind_rows(run_hyper("AEE"), run_hyper("PEE"))

enrich_AEE <- dplyr::bind_rows(run_hyper("AEE"))

enrich_PEE <- dplyr::bind_rows(run_hyper("PEE"))

sig_pathways <- enrich_all %>%
  dplyr::filter(p < 0.05, k >= 3) %>%
  dplyr::pull(Pathway) %>%
  unique()

enrich_all_filter <- enrich_all %>% subset(p < 0.05 & k >= 3)

enrich_all_filter <- left_join(enrich_all_filter, kegg_pathway_hierarchy, by = c('Pathway'='pathway_id'))

enrich_AEE_filter <- enrich_AEE %>% subset(p < 0.05 & k >= 3)

enrich_AEE_filter <- left_join(enrich_AEE_filter, kegg_pathway_hierarchy, by = c('Pathway'='pathway_id'))

enrich_PEE_filter <- enrich_PEE %>% subset(p < 0.05 & k >= 3)

enrich_PEE_filter <- left_join(enrich_PEE_filter, kegg_pathway_hierarchy, by = c('Pathway'='pathway_id'))

write.xlsx(enrich_all_filter, './微生物/KO_AEE_PEE/pathways enriched by differential ko.xlsx')


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ggforce)
  library(scales)
})

# =========================
# 0) parameters
# =========================
P_CUT <- 0.05
LABEL_LEVELS <- c("L1")   # 改成 c("L1","L2") 就会给 L2 也加标签
R0 <- 0
R1 <- 1.2
R2 <- 2.2
R3 <- 3.2

# =========================
# 2) KO -> Pathway -> Hierarchy (mapXXXXX only)
#    kegg_pathway_hierarchy: level1 level2 pathway_id name
# =========================
ko2hier <- ko_pathway %>%
  mutate(
    KO = sub("^ko:", "", KO),
    Pathway = as.character(Pathway)
  ) %>%
  filter(str_detect(Pathway, "^map\\d{5}$")) %>%
  inner_join(sig_ko, by = "KO") %>%
  filter(Pathway %in% sig_pathways) %>%   # ✅ 新增：只画显著富集的 pathway
  inner_join(
    kegg_pathway_hierarchy %>%
      transmute(Pathway = pathway_id, level1, level2, level3 = name),
    by = "Pathway"
  ) %>%
  mutate(
    level1 = ifelse(is.na(level1) | level1 == "", "Unclassified", level1),
    level2 = ifelse(is.na(level2) | level2 == "", "Unclassified", level2),
    level3 = ifelse(is.na(level3) | level3 == "", "Unclassified", level3)
  ) %>%
  distinct(KO, dir, Pathway, level1, level2, level3)

# ko2hier <- ko2hier %>% subset(!level1 %in% c('Human Diseases','Organismal Systems') )

# =========================
# 3) Count differential KO numbers at each node (L1/L2/L3) + add Center (L0)
# =========================
count_level3 <- ko2hier %>%
  group_by(level1, level2, level3, dir) %>%
  summarise(n = n_distinct(KO), .groups = "drop") %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(
    total = AEE + PEE, level = "L3",
    vid = paste0("L3__", level1, "__", level2, "__", level3),
    label = level3
  )

count_level2 <- ko2hier %>%
  group_by(level1, level2, dir) %>%
  summarise(n = n_distinct(KO), .groups = "drop") %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(
    total = AEE + PEE, level = "L2",
    vid = paste0("L2__", level1, "__", level2),
    label = level2
  )

count_level1 <- ko2hier %>%
  group_by(level1, dir) %>%
  summarise(n = n_distinct(KO), .groups = "drop") %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(
    total = AEE + PEE, level = "L1",
    vid = paste0("L1__", level1),
    label = level1
  )

count_center <- sig_ko %>%
  dplyr::count(dir, name = "n") %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(
    total = AEE + PEE, level = "L0",
    level1 = NA_character_, level2 = NA_character_, level3 = NA_character_,
    vid = "L0__KEGG",
    label = "KEGG"
  ) %>%
  dplyr::select(level1, level2, level3, AEE, PEE, total, level, vid, label)

nodes_df <- bind_rows(
  count_center %>%
    dplyr::select(level1, level2, level3, AEE, PEE, total, level, vid, label),
  
  count_level1 %>%
    dplyr::mutate(
      level2 = NA_character_,
      level3 = NA_character_
    ) %>%
    dplyr::select(level1, level2, level3, AEE, PEE, total, level, vid, label),
  
  count_level2 %>%
    dplyr::mutate(
      level3 = NA_character_
    ) %>%
    dplyr::select(level1, level2, level3, AEE, PEE, total, level, vid, label),
  
  count_level3 %>%
    dplyr::select(level1, level2, level3, AEE, PEE, total, level, vid, label)
) %>%
  dplyr::filter(total > 0) %>%
  dplyr::distinct(vid, .keep_all = TRUE)

# =========================
# 4) Build edges (L0->L1, L1->L2, L2->L3)
# =========================
edges_0_1 <- nodes_df %>%
  filter(level == "L1") %>%
  transmute(from = "L0__KEGG", to = vid)

edges_1_2 <- ko2hier %>%
  distinct(level1, level2) %>%
  transmute(
    from = paste0("L1__", level1),
    to   = paste0("L2__", level1, "__", level2)
  )

edges_2_3 <- ko2hier %>%
  distinct(level1, level2, level3) %>%
  transmute(
    from = paste0("L2__", level1, "__", level2),
    to   = paste0("L3__", level1, "__", level2, "__", level3)
  )

edges_df <- bind_rows(edges_0_1, edges_1_2, edges_2_3) %>%
  distinct() %>%
  filter(from %in% nodes_df$vid, to %in% nodes_df$vid)

# =========================
# 5) Custom radial layout (center -> L1 ring -> L2 ring -> L3 ring)
#    Angle allocation proportional to "total" so big categories get more space
# =========================

# helper: allocate [start,end] angle per child within a parent span, proportional to weight
allocate_span <- function(df, id_col, w_col, start, end) {
  # 防御：start/end 为空或 NA，就给整圈
  if (length(start) == 0 || length(end) == 0 || is.na(start) || is.na(end)) {
    start <- 0
    end <- 2*pi
  }
  if (end <= start) end <- start + 2*pi
  
  df <- df %>% arrange(desc(.data[[w_col]]))
  w <- df[[w_col]]
  w <- ifelse(is.na(w) | w <= 0, 1, w)
  
  cum <- cumsum(w) / sum(w)
  a_end <- start + (end - start) * cum
  a_start <- c(start, head(a_end, -1))
  
  df %>% mutate(theta_start = a_start, theta_end = a_end, theta = (a_start + a_end) / 2)
}

# 建议统一清洗 level1/level2/level3，避免不可见空格导致 join 不上
layout_df <- nodes_df %>%
  mutate(
    level1 = ifelse(is.na(level1), NA_character_, str_squish(level1)),
    level2 = ifelse(is.na(level2), NA_character_, str_squish(level2)),
    level3 = ifelse(is.na(level3), NA_character_, str_squish(level3))
  ) %>%
  mutate(r_level = case_when(level=="L0" ~ R0, level=="L1" ~ R1, level=="L2" ~ R2, level=="L3" ~ R3))

# L1: full circle（一定要保留 theta_start/theta_end）
l1 <- layout_df %>% filter(level=="L1")
l1 <- allocate_span(l1, "vid", "total", 0, 2*pi) %>%
  mutate(x = R1 * cos(theta), y = R1 * sin(theta)) %>%
  select(vid, level1, theta_start, theta_end, theta, x, y, everything())

# L2: within each L1 span
l2 <- layout_df %>% filter(level=="L2") %>%
  left_join(l1 %>% distinct(level1, theta_start, theta_end), by="level1") %>%
  mutate(
    theta_start = ifelse(is.na(theta_start), 0, theta_start),
    theta_end   = ifelse(is.na(theta_end), 2*pi, theta_end)
  ) %>%
  group_by(level1) %>%
  group_modify(~ allocate_span(.x, "vid", "total", unique(.x$theta_start)[1], unique(.x$theta_end)[1])) %>%
  ungroup() %>%
  mutate(x = R2 * cos(theta), y = R2 * sin(theta))

# L3: within each L2 span
l3 <- layout_df %>% filter(level=="L3") %>%
  left_join(l2 %>% select(level1, level2, theta_start, theta_end), by=c("level1","level2")) %>%
  group_by(level1, level2) %>%
  group_modify(~ allocate_span(.x, "vid", "total", unique(.x$theta_start), unique(.x$theta_end))) %>%
  ungroup() %>%
  mutate(x = R3 * cos(theta), y = R3 * sin(theta))

# merge back
layout_df <- layout_df %>%
  dplyr::select(-dplyr::any_of(c("x","y","theta","theta_start","theta_end"))) %>%
  dplyr::left_join(
    dplyr::bind_rows(
      layout_df %>%
        dplyr::filter(level == "L0") %>%
        dplyr::mutate(theta = 0, theta_start = 0, theta_end = 2*pi, x = 0, y = 0) %>%
        dplyr::select(vid, theta, theta_start, theta_end, x, y),
      
      l1 %>% dplyr::select(vid, theta, theta_start, theta_end, x, y),
      l2 %>% dplyr::select(vid, theta, theta_start, theta_end, x, y),
      l3 %>% dplyr::select(vid, theta, theta_start, theta_end, x, y)
    ) %>% dplyr::distinct(vid, .keep_all = TRUE),
    by = "vid"
  )

# =========================
# 6) Build graph + attach manual coords to layout
# =========================
vertices_df <- nodes_df %>%
  transmute(name = vid, level, label, AEE, PEE, total) %>%
  distinct(name, .keep_all = TRUE) %>%
  select(name, everything())

# =========================
# Tree-only plot (no pies)
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
})

# 1) build graph
# vertices_df 必须第一列叫 name（唯一ID），且包含 level/label/total/AEE/PEE
g <- igraph::graph_from_data_frame(edges_df, vertices = vertices_df, directed = TRUE)

# 2) tree plot
p_tree <- ggraph(g, layout = "dendrogram",circular = TRUE) +
  geom_edge_diagonal(colour = "black", linewidth = 0.6) +
  geom_node_point(
    aes(size = total),
    colour = "grey20",
    alpha = 1
  ) +
  # L3 全标注（你要求的三级通路全部加名字）
  geom_node_text(
    data = function(x) subset(x, level == "L3"),
    aes(label = label),
    size = 2.6,
    colour = "black",
    hjust = 0,
    nudge_x = 0.02
  ) +
  # 可选：L1/L2 也加（建议小一点，避免挤）
  geom_node_text(
    data = function(x) subset(x, level %in% c("L1","L2")),
    aes(label = label),
    size = 2.4,
    colour = "black",
    vjust = -0.8
  ) +
  scale_size(range = c(6,20)) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
  ) +
  labs(size = "Enriched KO count")

p_tree

ggsave('./微生物/KO_AEE_PEE_CLR_LIMMA/AEE_PEE_diff_KO_pathway.pdf',p_tree, width = 13, height = 6)




# =========================
# One pie per node (export)
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

OUT_PIE_DIR <- "./微生物/KO_AEE_PEE_CLR_LIMMA/pies_by_node"
dir.create(OUT_PIE_DIR, showWarnings = FALSE, recursive = TRUE)

# 只给 total>0 的节点输出饼图（建议：只输出 L3 或者你关心的层级）
pie_nodes <- nodes_df %>%
  filter(total > 0) %>%
  # 你要“每个节点一张图”，如果太多建议先只做 L3：
  # filter(level == "L3") %>%
  select(vid, label, level, AEE, PEE, total)

make_one_pie <- function(label, aee, pee, out_png) {
  df <- tibble(
    group = c("AEE","PEE"),
    n = c(aee, pee)
  ) #%>% filter(n > 0)
  
  # 如果某个节点只有单边（aee=0 或 pee=0），也能画
  p <- ggplot(df, aes(x = 1, y = n, fill = group)) +
    geom_col(width = 1, colour = NA) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(AEE = "#2e528f", PEE = "#b41d23")) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA)
    )
  
  ggsave(out_png, p, width = 2, height = 2, dpi = 300, bg = "transparent")
}

pie_nodes$label <- gsub(": ", "", pie_nodes$label)

for (i in seq_len(nrow(pie_nodes))) {
  vid_i   <- pie_nodes$vid[i]
  lab_i   <- pie_nodes$label[i]
  lev_i   <- pie_nodes$level[i]
  aee_i   <- pie_nodes$AEE[i]
  pee_i   <- pie_nodes$PEE[i]
  
  # 文件名要安全（去掉特殊字符）
  safe_lab <- lab_i %>%
    gsub("[^A-Za-z0-9_\\-]+", "_", .) %>%  # 非法字符 → _
    gsub("_+", "_", .) %>%                # 多个 _ 合并
    gsub("^_|_$", "", .) %>% 
    gsub(": ", "", .)
  
  out_png <- file.path(
    OUT_PIE_DIR,
    sprintf("%s__%s__%s.pdf", lev_i, safe_lab, vid_i)
  )
  
  make_one_pie(lab_i, aee_i, pee_i, out_png)
}

message("Done. Pies saved to: ", OUT_PIE_DIR)

i=29
vid_i   <- pie_nodes$vid[i]
lab_i   <- pie_nodes$label[i]
lev_i   <- pie_nodes$level[i]
aee_i   <- pie_nodes$AEE[i]
pee_i   <- pie_nodes$PEE[i]

# 文件名要安全（去掉特殊字符）
safe_lab <- lab_i %>%
  gsub("[^A-Za-z0-9_\\-]+", "_", .) %>%  # 非法字符 → _
  gsub("_+", "_", .) %>%                # 多个 _ 合并
  gsub("^_|_$", "", .)                  # 去掉首尾 _
out_png <- file.path(OUT_PIE_DIR, sprintf("%s__%s__%s.pdf", lev_i, safe_lab, vid_i))

df <- tibble(
  group = c("AEE","PEE"),
  n = c(aee_i, pee_i)
) 

# 如果某个节点只有单边（aee=0 或 pee=0），也能画
ggplot(df, aes(x = 1, y = n, fill = group)) +
  geom_col(width = 1, colour = NA) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(AEE = "#2e528f", PEE = "#b41d23")) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

ggsave(out_png,width = 2, height = 2, dpi = 300, bg = "transparent")

#### AEE PEE 分开看
# =========================
# AEE Hypergeometric enrichment on KEGG pathways
# =========================
P_CUT <- 0.05
ENRICH_FDR <- 0.1
MIN_K <- 3

df_merge$KO <- sub('^ko:','',df_merge$KO)

# universe
universe_ko <- df_merge$KO

# nominal diff KO
sig_ko_AEE <- df_merge %>%
  dplyr::mutate(
    KO  = sub("^ko:", "", KO),
    dir = dplyr::case_when(
      beta_AP_CLR*beta_group_CLR > 0 & beta_AP_CLR > 0 & (pvalue_AP < P_CUT| pvalue_group< P_CUT)  ~ "AEE",
      beta_AP_CLR*beta_group_CLR > 0 & beta_AP_CLR < 0 & (pvalue_AP < P_CUT| pvalue_group< P_CUT)  ~ "PEE",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(dir), pvalue_AP < P_CUT | pvalue_group < P_CUT, dir == 'AEE') %>%
  dplyr::select(KO, dir) %>%
  dplyr::distinct()

# KO -> Pathway (map only) within universe
ko2path_all <- ko_pathway %>%
  dplyr::mutate(
    KO = sub("^ko:", "", KO),
    Pathway = as.character(Pathway)
  ) %>%
  dplyr::filter(stringr::str_detect(Pathway, "^map\\d{5}$")) %>%
  dplyr::filter(KO %in% universe_ko) %>%
  dplyr::distinct(KO, Pathway)

universe_mappable <- unique(ko2path_all$KO)
N <- length(universe_mappable)

path_bg <- ko2path_all %>%
  dplyr::count(Pathway, name = "M")

run_hyper <- function(dir_label) {
  query_ko <- sig_ko_AEE %>%
    dplyr::filter(dir == dir_label) %>%
    dplyr::pull(KO) %>%
    unique()
  
  query_ko <- intersect(query_ko, universe_mappable)
  n <- length(query_ko)
  
  hit <- ko2path_all %>%
    dplyr::filter(KO %in% query_ko) %>%
    dplyr::count(Pathway, name = "k")
  
  path_bg %>%
    dplyr::left_join(hit, by = "Pathway") %>%
    dplyr::mutate(
      k = tidyr::replace_na(k, 0L),
      n = n,
      N = N,
      p = stats::phyper(q = k - 1, m = M, n = N - M, k = n, lower.tail = FALSE),
      dir = dir_label
    ) %>%
    dplyr::mutate(padj = p.adjust(p, method = "BH")) %>%
    dplyr::arrange(p)
}

enrich_all <- dplyr::bind_rows(run_hyper("AEE"), run_hyper("PEE"))

sig_pathways <- enrich_all %>%
  dplyr::filter(p < 0.05, k >= 3) %>%
  dplyr::pull(Pathway) %>%
  unique()

enrich_all_filter <- enrich_all %>% subset(p < 0.05 & k >= 3)

enrich_all_filter <- left_join(enrich_all_filter, kegg_pathway_hierarchy, by = c('Pathway'='pathway_id'))

write.xlsx(enrich_all_filter, './微生物/KO_AEE_PEE/pathways enriched by differential ko.xlsx')



