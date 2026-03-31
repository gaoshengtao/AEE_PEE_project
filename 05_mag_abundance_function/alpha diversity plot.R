# =========================
# LEfSe input from CoverM + rumenotype assignments
# + Alpha diversity + group test + plot
# + NEW: LDA top10 MAGs (non-Others Species) abundance ~ AP regression per rumenotype
# =========================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(lvplot)
  library(ggsignif)
  library(colorspace)
  library(grid)
  
  # you used openxlsx
  library(openxlsx)
  library(patchwork)   # 拼小图用（没有也行，可注释掉，改成逐个保存）
})

# ---- input paths ----
coverm_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/repMAGs_coverm.tsv"
assign_xlsx <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx"

# NEW: differential MAGs excel
diffmag_xlsx <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 4-AEE PEE differential MAGs.xlsx"

# ---- outputs ----
out_alpha_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/alpha_diversity_by_rumenotype.tsv"
out_alpha_p_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/alpha_diversity_group_tests.tsv"
out_alpha_pdf <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/alpha_diversity_lvbox.pdf"

# NEW outputs
out_map_xlsx <- "./微生物/VFA-metagenome sampleID mapping.xlsx"
out_lm_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/mag_ap_lm_results.tsv"
out_lm_pdf_A <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/figures/MAG_AP_regression_top10_Acetate.pdf"
out_lm_pdf_P <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/figures/MAG_AP_regression_top10_Propionate.pdf"

# ---- parameters ----
abund_suffix <- "TPM"  # or "Relative Abundance (%)"
GROUP_KEEP <- c("Acetate_type", "Propionate_type")
TEST_METHOD <- "wilcox"

COLORS_BASE <- c(
  "Acetate_type"    = "#2F5597",
  "Propionate_type" = "#C00000"
)

TOP_N_LDA <- 10
SPECIES_OTHERS <- c("Others", "Other", "others", "other", "NA", "", NA)

# =========================
# 1) Read group assignments
# =========================
meta <- readxl::read_excel(assign_xlsx, sheet = 1) %>% as.data.frame()

VFA <- readxl::read_excel(assign_xlsx, sheet = 2) %>% as.data.frame()

meta <- left_join(meta, VFA, by='SampleID')

# 这里假设表里已经有 SampleID / Rumenotype / AP 列
meta <- meta %>%
  mutate(
    SampleID = as.character(SampleID),
    digits = str_extract(SampleID, "\\d+"),
    Rumenotype = as.character(Rumenotype)
  ) %>%
  filter(Rumenotype %in% GROUP_KEEP) %>%
  mutate(Rumenotype = factor(Rumenotype, levels = GROUP_KEEP))

if (nrow(meta) == 0) stop("未找到属于两型的样本，请检查 Rumenotype 拼写。")

# =========================
# 2) Read CoverM table
# =========================
coverm <- readr::read_tsv(coverm_tsv, show_col_types = FALSE)
if (!"MAG_formal" %in% colnames(coverm)) stop("CoverM 表未找到 MAG_formal 列。")

# =========================
# 3) Choose abundance columns
# =========================
abund_cols <- colnames(coverm)[str_detect(colnames(coverm), fixed(abund_suffix, ignore_case = FALSE))]
if (length(abund_cols) == 0) stop(paste0("未找到包含后缀 '", abund_suffix, "' 的列。"))

sample_codes <- abund_cols %>%
  str_replace("\\.repMAGs\\.sorted.*$", "")
sample_digits <- str_extract(sample_codes, "\\d+")

map_coverm <- tibble(
  coverm_col = abund_cols,
  coverm_code = sample_codes,
  digits = sample_digits
)

# =========================
# 4) Match samples by digits
# =========================
matched <- meta %>% inner_join(map_coverm, by = "digits")
if (nrow(matched) == 0) stop("meta 与 coverm 未匹配到样本（digits 不一致）。")

matched <- matched %>%
  mutate(meta_order = match(SampleID, meta$SampleID)) %>%
  arrange(meta_order)

# openxlsx::write.xlsx(matched, out_map_xlsx)

final_sample_ids  <- matched$SampleID
final_class       <- matched$Rumenotype
final_coverm_cols <- matched$coverm_col

# =========================
# 5) Build feature table (feat) + LEfSe input  ✅修复 feat 未定义 bug
# =========================
feat <- coverm %>%
  dplyr::select(MAG_formal, all_of(final_coverm_cols)) %>%
  mutate(MAG_formal = as.character(MAG_formal)) %>%
  filter(!is.na(MAG_formal), MAG_formal != "unmapped")   # 更稳健：直接去掉 unmapped 行

# numeric coerce
for (cc in final_coverm_cols) feat[[cc]] <- suppressWarnings(as.numeric(feat[[cc]]))
feat[is.na(feat)] <- 0

# =========================
# 6) Alpha diversity (per sample)
# =========================
X <- as.matrix(feat[, final_coverm_cols, drop = FALSE])
rownames(X) <- feat$MAG_formal

X_samp <- t(X)                         # rows=sample, cols=genome
rownames(X_samp) <- final_sample_ids

Observed <- rowSums(X_samp > 0)
Shannon  <- vegan::diversity(X_samp, index = "shannon")
Simpson  <- vegan::diversity(X_samp, index = "simpson")

# Chao1：estimateR 用 counts 更合理，你用了 TPM，这里按你的思路 round()
estR <- vegan::estimateR(t(round(X_samp)))
Chao1 <- as.numeric(estR["S.chao1", ])
names(Chao1) <- colnames(t(round(X_samp)))

alpha_diversity <- tibble(
  SampleID = final_sample_ids,
  Rumenotype = as.character(final_class),
  Observed = as.numeric(Observed),
  Chao1 = as.numeric(Chao1[final_sample_ids]),
  Shannon = as.numeric(Shannon),
  Simpson = as.numeric(Simpson)
) %>%
  left_join(meta, by = c("SampleID", "Rumenotype")) %>%
  mutate(Rumenotype = factor(Rumenotype, levels = GROUP_KEEP))

readr::write_tsv(alpha_diversity, out_alpha_tsv)
message("✅ alpha 多样性表已输出：", out_alpha_tsv)

# =========================
# 7) Group test
# =========================
run_test <- function(df, metric, group_col = "Rumenotype", method = c("wilcox", "t.test")) {
  method <- match.arg(method)
  df <- df %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[metric]]))
  if (n_distinct(df[[group_col]]) != 2) return(NA_real_)
  if (method == "wilcox") suppressWarnings(wilcox.test(df[[metric]] ~ df[[group_col]])$p.value)
  else suppressWarnings(t.test(df[[metric]] ~ df[[group_col]])$p.value)
}

pvals <- tibble(
  metric = c("Observed", "Chao1", "Shannon", "Simpson"),
  p_value = c(
    run_test(alpha_diversity, "Observed", method = TEST_METHOD),
    run_test(alpha_diversity, "Chao1", method = TEST_METHOD),
    run_test(alpha_diversity, "Shannon", method = TEST_METHOD),
    run_test(alpha_diversity, "Simpson", method = TEST_METHOD)
  ),
  method = TEST_METHOD
)

readr::write_tsv(pvals, out_alpha_p_tsv)
message("✅ 组间差异检验结果已输出：", out_alpha_p_tsv)

# =========================
# 8) Alpha plot (your style)
# =========================
alpha_long <- alpha_diversity %>%
  pivot_longer(cols = c(Observed, Chao1, Shannon, Simpson),
               names_to = "metric", values_to = "value")

p_all <- ggplot(alpha_long, aes(x = Rumenotype, y = value, fill = Rumenotype, colour = Rumenotype)) +
  lvplot::geom_lv(k = 5, outlier.shape = NA) +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = "#00000000", size = 1.2) +
  geom_jitter(shape = 21, size = 2.8, fill = "white", width = 0.1, height = 0) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(2, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  scale_color_manual(values = setNames(darken(COLORS_BASE, 0.5), names(COLORS_BASE))) +
  scale_fill_manual(values  = setNames(darken(COLORS_BASE, 0.0), names(COLORS_BASE))) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  geom_signif(comparisons = list(GROUP_KEEP), color = "black", tip_length = 0) +
  labs(x = NULL, y = NULL)

ggsave(out_alpha_pdf, p_all, width = 10, height = 4)
message("✅ alpha 多样性图已输出：", out_alpha_pdf)

# =========================
# 9) NEW: read differential MAGs, pick top10 per rumenotype, regression abundance ~ AP
# =========================
diffmag <- readxl::read_excel(diffmag_xlsx, sheet = 1) %>% as.data.frame()

need_cols <- c("MAG_formal","Rumenotype","LDA_score","pvalue","Species")
miss <- setdiff(need_cols, colnames(diffmag))
if (length(miss) > 0) stop("Supplemental file 4 缺少列：", paste(miss, collapse = ", "))

diffmag <- diffmag %>%
  mutate(
    MAG_formal = as.character(MAG_formal),
    Rumenotype = as.character(Rumenotype),
    LDA_score = suppressWarnings(as.numeric(LDA_score)),
    pvalue = suppressWarnings(as.numeric(pvalue)),
    Species = as.character(Species)
  ) %>%
  filter(Rumenotype %in% c("Acetate_type","Propionate_type")) %>%
  filter(!(Species %in% SPECIES_OTHERS))

# 每个组取 LDA top10
# 如果你的 LDA 可能有负值，想按 |LDA| 排序可改成 arrange(desc(abs(LDA_score)))
top_mag <- diffmag %>%
  group_by(Rumenotype) %>%
  arrange(desc(LDA_score)) %>%
  slice_head(n = TOP_N_LDA) %>%
  ungroup()

if (nrow(top_mag) == 0) stop("在 diff MAG 表里没有筛到 top MAG（检查 Rumenotype / Species 列值）。")

# --- 构建 abundance long table：每行=SampleID+MAG ---
abund_long <- feat %>%
  select(MAG_formal, all_of(final_coverm_cols)) %>%
  pivot_longer(cols = all_of(final_coverm_cols), names_to = "coverm_col", values_to = "abundance") %>%
  left_join(matched, by = "coverm_col") %>%
  rename(MAG_formal = MAG_formal) %>%
  filter(!is.na(SampleID), !is.na(Rumenotype))

# 只保留 top MAG
abund_top <- abund_long %>%
  inner_join(top_mag %>% select(MAG_formal, Rumenotype, LDA_score, pvalue, Species),
             by = c("MAG_formal", "Rumenotype"))

Varable = 

# --- 逐 MAG 做线性回归：abundance ~ AP（在各自组内）---
fit_one <- function(df){
  df <- df 
  if (nrow(df) < 3) return(tibble(n = nrow(df), r2 = NA_real_, p = NA_real_))
  m <- lm(abundance ~ Acetic_acid, data = df)
  sm <- summary(m)
  tibble(
    n = nrow(df),
    r2 = unname(sm$r.squared),
    p  = unname(coef(sm)[2,4])   # AP 的 p 值
  )
}

lm_res <- abund_top %>%
  group_by(Rumenotype, MAG_formal, Species, LDA_score, pvalue) %>%
  group_modify(~fit_one(.x)) %>%
  ungroup() %>%
  arrange(Rumenotype, desc(LDA_score))
# 
# readr::write_tsv(lm_res, out_lm_tsv)
# message("✅ MAG ~ AP 回归结果已输出：", out_lm_tsv)

# =========================
# 10) Plot small multiples with R2 and P on each panel
# =========================
# 先把回归结果合并回点数据，方便 annotate
abund_top2 <- abund_top %>%
  left_join(lm_res, by = c("Rumenotype","MAG_formal","Species","LDA_score","pvalue"))


make_panel <- function(df_one){
  # df_one: 单个 MAG 的数据
  lab <- unique(df_one$MAG_formal)
  sp  <- unique(df_one$Species)
  r2  <- unique(df_one$r2)
  pv  <- unique(df_one$p)
  
  ann <- paste0("R² = ", ifelse(is.na(r2), "NA", sprintf("%.2f", r2)),
                "\nP = ",  ifelse(is.na(pv), "NA", format.pval(pv, digits = 2, eps = 1e-3)))
  
  ggplot(df_one, aes(x = Acetic_acid, y = abundance)) +
    geom_point(shape = 21, size = 2.5, fill = "white", colour = "black") +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 10)
    ) +
    labs(
      title = paste0(lab, " | ", sp),
      x = "Acetic_acid",
      y = abund_suffix
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.1,
             label = ann, size = 3)
}

plot_group_pdf <- function(group_name, out_pdf){
  df_g <- abund_top2 %>% filter(Rumenotype == group_name)
  
  mags <- unique(df_g$MAG_formal)
  if (length(mags) == 0) {
    message("⚠️ 组 ", group_name, " 没有可画的 MAG（可能都被过滤或未匹配到丰度）。")
    return(invisible(NULL))
  }
  
  plist <- lapply(mags, function(mg){
    make_panel(df_g %>% filter(MAG_formal == mg))
  })
  
  # 拼图：2行5列（top10）
  wrap <- patchwork::wrap_plots(plist, ncol = 5)
  ggsave(out_pdf, wrap, width = 16, height = 6)
  message("✅ ", group_name, " 的 top MAG 回归小图已输出：", out_pdf)
}

plot_group_pdf("Acetate_type", out_lm_pdf_A)
plot_group_pdf("Propionate_type", out_lm_pdf_P)


# Propionic_acid
make_panel <- function(df_one){
  # df_one: 单个 MAG 的数据
  lab <- unique(df_one$MAG_formal)
  sp  <- unique(df_one$Species)
  r2  <- unique(df_one$r2)
  pv  <- unique(df_one$p)
  
  ann <- paste0("R² = ", ifelse(is.na(r2), "NA", sprintf("%.2f", r2)),
                "\nP = ",  ifelse(is.na(pv), "NA", format.pval(pv, digits = 2, eps = 1e-3)))
  
  ggplot(df_one, aes(x = Propionic_acid, y = abundance)) +
    geom_point(shape = 21, size = 2.5, fill = "white", colour = "black") +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 10)
    ) +
    labs(
      title = paste0(lab, " | ", sp),
      x = "Propionic_acid",
      y = abund_suffix
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.1,
             label = ann, size = 3)
}

plot_group_pdf <- function(group_name, out_pdf){
  df_g <- abund_top2 %>% filter(Rumenotype == group_name)
  
  mags <- unique(df_g$MAG_formal)
  if (length(mags) == 0) {
    message("⚠️ 组 ", group_name, " 没有可画的 MAG（可能都被过滤或未匹配到丰度）。")
    return(invisible(NULL))
  }
  
  plist <- lapply(mags, function(mg){
    make_panel(df_g %>% filter(MAG_formal == mg))
  })
  
  # 拼图：2行5列（top10）
  wrap <- patchwork::wrap_plots(plist, ncol = 5)
  ggsave(out_pdf, wrap, width = 16, height = 6)
  message("✅ ", group_name, " 的 top MAG 回归小图已输出：", out_pdf)
}

plot_group_pdf("Acetate_type", out_lm_pdf_A)
plot_group_pdf("Propionate_type", out_lm_pdf_P)




