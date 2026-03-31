# =========================
# mo-based functional analysis: PCoA + adonis2 + DESeq2 + Volcano
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(DESeq2)
  library(ggrepel)
})

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(ggplot2)
})


# --------- 0) Inputs (EDIT) ----------
mo_tsv   <- "./微生物/genecatelog/06_aggregate_tpm_by_eggnog/TPM_by_KEGG_Module.tsv"   # mo x Sample matrix (tsv)
meta_xlsx <- "./FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx"            # Sample metadata with SampleID + Group (AEE/PEE)
outdir   <- "./微生物/mo_AEE_PEE"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

module_name <- read.table("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/module name.txt", header = T, sep = '\t')

kegg_module_hierarchy <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/kegg_module_hierarchy_L1L2L3.tsv")

infer_module <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/genecatelog/06_infer_kegg_names_final/modules_bestname.tsv")

infer_module <- left_join(infer_module, kegg_module_hierarchy, by = c('Best_KEGG_Module' = 'module_id'))



# metadata columns
sample_col <- "SampleID"
group_col  <- "Rumenotype"   # levels: Acetate_type / Propionate_type
ap_col     <- "AP"

# filtering & CLR params
min_prevalence <- 0.10  # keep mo present in >=10% samples (TPM > 0)
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
# 2) Read mo TPM matrix
# ------------------------------------------------------------
mo <- fread(mo_tsv)
names(mo)[1] <- "mo"

sample_cols <- intersect(names(mo), meta[[sample_col]])
stopifnot(length(sample_cols) >= 4)

mo_mat <- mo %>%
  select(mo, all_of(sample_cols)) %>%
  as.data.frame()

rownames(mo_mat) <- mo_mat$mo
mo_mat$mo <- NULL

# align meta to matrix columns
meta2 <- meta %>%
  filter(.data[[sample_col]] %in% colnames(mo_mat)) %>%
  arrange(match(.data[[sample_col]], colnames(mo_mat)))

mo_mat <- mo_mat[, meta2[[sample_col]], drop = FALSE]

stopifnot(all(colnames(mo_mat) == meta2[[sample_col]]))

# ------------------------------------------------------------
# 3) Filter low-prevalence mos (based on TPM > 0)
# ------------------------------------------------------------
prev <- rowMeans(mo_mat > 0, na.rm = TRUE)
keep <- prev >= min_prevalence
mo_f <- mo_mat[keep, , drop = FALSE]

message(sprintf("[filter] mos: %d (raw) -> %d (prevalence >= %.0f%%)",
                nrow(mo_mat), nrow(mo_f), 100 * min_prevalence))

# ------------------------------------------------------------
# 4) TPM -> CLR
# ------------------------------------------------------------
clr_mat <- clr_transform(as.matrix(mo_f), pseudo = pseudo_count)

# Remove rows with zero variance (can cause NA coef)
mo_var <- apply(clr_mat, 1, var, na.rm = TRUE)
clr_mat <- clr_mat[mo_var > 0, , drop = FALSE]

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
  tibble::rownames_to_column("mo") %>%
  as_tibble() %>%
  transmute(
    mo,
    beta_AP_CLR = logFC,
    pvalue_AP   = P.Value,
    padj_AP     = adj.P.Val
  ) %>%
  arrange(padj_AP, pvalue_AP)

res_ap$color <- ifelse(res_ap$beta_AP_CLR > 0 & res_ap$pvalue_AP < 0.05, 'AP_pos',
                       ifelse(res_ap$beta_AP_CLR < 0 & res_ap$pvalue_AP < 0.05, 'AP_neg', 'N.S.'))

safe_write(res_ap, file.path(outdir, "mo_betaAP_CLR_limma.tsv"))

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

ggsave(file.path(outdir, "mo_betaAP_CLR_vol_scatter.pdf"), p_ap, width = 4, height = 4)

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
  tibble::rownames_to_column("mo") %>%
  as_tibble()

res_g <- tmp_g %>%
  transmute(
    mo,
    beta_group_CLR = -logFC,     # Acetate - Propionate
    pvalue_group   = P.Value,
    padj_group     = adj.P.Val
  ) %>%
  arrange(padj_group, pvalue_group)

res_g$color <- ifelse(res_g$beta_group_CLR > 0 & res_g$pvalue_group < 0.05, 'AP_pos',
                      ifelse(res_g$beta_group_CLR < 0 & res_g$pvalue_group < 0.05, 'AP_neg', 'N.S.'))

safe_write(res_g, file.path(outdir, "mo_betaGroup_CLR_AminusP_limma.tsv"))

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

ggsave(file.path(outdir, "mo_betagroup_CLR_vol_scatter.pdf"), p_g, width = 4, height = 4)

table(res_g$color)

# ------------------------------------------------------------
# 7) Consistency: correlation + quadrant plot (same-sign fraction)
# ------------------------------------------------------------
df_merge <- res_ap %>%
  inner_join(res_g, by = "mo") %>%
  mutate(
    consistent = (beta_AP_CLR * beta_group_CLR) > 0,
    cons_lab = ifelse(consistent, "Consistent", "Inconsistent")
  )

df_merge$cons_lab <- ifelse(df_merge$color.x == 'N.S.' & df_merge$color.y == 'N.S.', 'N.S.', df_merge$cons_lab)

safe_write(df_merge, file.path(outdir, "mo_betaAP_vs_betaGroup_CLR_merged.tsv"))

cons_summary <- df_merge %>%
  summarise(
    n = n(),
    prop_consistent = mean(consistent, na.rm = TRUE),
    spearman = cor(beta_AP_CLR, beta_group_CLR, method = "spearman", use = "pairwise.complete.obs"),
    pearson  = cor(beta_AP_CLR, beta_group_CLR, method = "pearson",  use = "pairwise.complete.obs")
  )
safe_write(cons_summary, file.path(outdir, "consistency_summary.tsv"))

COL_CONS <- c("Consistent" = "#D62728", "Inconsistent" = "#1F77B4", "N.S."='#efefef')

df_merge <- left_join(df_merge, infer_module, by= c('mo'='EggNOG_Module'))

df_merge$quad_color <- 'N.S.'

df_merge$quad_color <- ifelse(df_merge$cons_lab == "Consistent", df_merge$level2,
                              ifelse(df_merge$cons_lab == "Inconsistent", 'Inconsistent', 'N.S.')
                              )
           
df_merge$quad_color[is.na(df_merge$quad_color)] <- 'Others'

module_colors <- c(
  # ---- Core metabolism modules ----
  "Carbohydrate metabolism"                    = "#4E79A7",  # muted blue
  "Energy metabolism"                          = "#59A14F",  # muted green
  "Lipid metabolism"                           = "#E15759",  # soft red
  "Amino acid metabolism"                      = "#76B7B2",  # teal
  "Nucleotide metabolism"                      = "#EDC948",  # muted yellow
  "Glycan metabolism"                          = "#B07AA1",  # purple
  "Metabolism of cofactors and vitamins"       = "#F28E2B",  # orange
  "Biosynthesis of terpenoids and polyketides" = "#9C755F",  # brown
  "Xenobiotics biodegradation"                 = "#8CD17D",  # light green
  
  # ---- Non-core / annotation states (grey gradient) ----
  "Gene set"                                   = "#C7ABCD",  # 深灰
  "Others"                                     = "#9E9E9E",  # 中灰
  "Inconsistent"                               = "#C7C7C7",  # 浅灰
  "N.S."                                       = "#EFEFEF"   # 极浅灰（你指定）
)

imposable <- c('Complete nitrification, comammox, ammonia => nitrite => nitrate','V-type ATPase, eukaryotes','Anoxygenic photosystem I',
               'N-glycan biosynthesis, complex type','Ascorbate degradation, ascorbate => D-xylulose-5P',
               'Castasterone biosynthesis, plants, campesterol => castasterone','GABA biosynthesis, eukaryotes, putrescine => GABA',
               'Neocarzinostatin naphthoate moiety biosynthesis','Cationic antimicrobial peptide (CAMP) resistance, protease PgtE',
               'C19-Steroid hormone biosynthesis, pregnenolone => testosterone => dihydrotestosterone')

df_merge[df_merge$Best_Name %in% imposable, ]$quad_color <- 'Others'
           
p_quad <- ggplot(df_merge, aes(x = beta_group_CLR, y = beta_AP_CLR, color = quad_color)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey35") +
  geom_point(alpha = 1, size = 3) +
  scale_color_manual(values = module_colors) +
  theme_bw(base_size = 11)+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm')
  )

p_quad

ggsave(file.path(outdir, "mo_betaAP_vs_betaGroup_CLR_quadrant.pdf"),
       p_quad, width = 4, height = 5.0)

df_merge[df_merge$beta_AP_CLR > 0 & ! df_merge$quad_color %in% c('Inconsistent','Others','N.S.'),]$Best_Name 

df_merge[df_merge$beta_AP_CLR < 0 & ! df_merge$quad_color %in% c('Inconsistent','Others','N.S.'),]$Best_Name 
# 

df_merge_circle <- df_merge %>% subset(! quad_color %in% c('Inconsistent','Others','N.S.'))


df_merge_circle$dir <- ifelse(df_merge_circle$beta_AP_CLR > 0, 'AEE', 'PEE')

table(df_merge_circle$dir, df_merge_circle$quad_color)

df_merge_circle$module_name

library(tidyverse)

# 你现有的 table
tab <- table(df_merge_circle$dir, df_merge_circle$quad_color)

df_bubble <- as.data.frame(tab) %>%
  dplyr::rename(
    Dir      = Var1,
    Category = Var2,
    Count    = Freq
  ) %>%
  filter(Count > 0) %>%        # 不画 0
  mutate(
    Dir = factor(Dir, levels = c("AEE", "PEE")),
    Category = factor(Category, levels = rev(unique(Category)))
  )

df_bubble <- df_bubble %>%
  group_by(Category) %>%
  mutate(total = sum(Count)) %>%
  ungroup() %>%
  mutate(Category = fct_reorder(Category, total))


p <- ggplot(df_bubble, aes(
  x = Dir,
  y = Category,
  size = Count,
  fill = Dir
)) +
  geom_point(
    shape = 21,
    color = "black",
    alpha = 0.85,
    stroke = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "AEE" = "#2e528f",
      "PEE" = "#b41d23"
    )
  ) +
  scale_size_continuous(
    range = c(3, 10)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Number of modules"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

p

ggsave(file.path(outdir, "mo_betaAP_vs_betaGroup_CLR_quadrant_circle.pdf"),
       p, width = 6, height = 5.0)



cons_order <- c(
  "Amino acid metabolism",
  "Carbohydrate metabolism",
  "Energy metabolism",
  "Glycan metabolism",
  "Gene set",
  "Biosynthesis of terpenoids and polyketides",
  "Lipid metabolism",
  "Metabolism of cofactors and vitamins",
  "Nucleotide metabolism",
  "Xenobiotics biodegradation"
)

df_merge_circle <- df_merge_circle %>%
  mutate(quad_color = factor(quad_color, levels = cons_order)) %>%
  arrange(quad_color, beta_AP_CLR)


df_merge_circle$Best_Name <- factor(df_merge_circle$Best_Name,levels = unique(df_merge_circle$Best_Name))

p <- ggplot(df_merge_circle, aes(
  x = beta_AP_CLR,
  y = Best_Name,
  color = quad_color
)) +
  # 棒子
  geom_segment(
    aes(x = 0, xend = beta_AP_CLR, yend = Best_Name),
    linewidth = 1.2,
    alpha = 0.8
  ) +
  # 棒头
  geom_point(
    size = 3.5
  ) +
  scale_color_manual(
    values = module_colors
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(5),
    labels = abs
  ) +
  labs(
    x = "Number of enriched modules",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'none'
  )

p


ggsave(file.path(outdir, "mo_betaAP_vs_betaGroup_CLR_quadrant_circle_bbt.pdf"),
       p, width = 12, height = 5.0)


























