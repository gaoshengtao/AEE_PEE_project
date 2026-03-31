# ============================================================
# Barplot: LEfSe MAGs (LDA) by Species (join by MAG_formal)
# ============================================================

pkgs <- c("readxl", "dplyr", "ggplot2", "stringr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
lapply(pkgs, library, character.only = TRUE)

# ---- paths ----
tax_path   <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 6-master_repMAG_taxonomy_by_rank.xlsx"
lefse_path <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/lefse.LDA_repMAGs_tpm.xlsx"

# ---- colors ----
type_cols <- c(
  "Acetate_type"    = "#2F5597",
  "Propionate_type" = "#C00000"
)

# ---- parameters (你可按需要改) ----
alpha_p <- 0.05          # 显著性阈值（如果你不想按P筛，就改成 1）
lda_cut <- 0             # 可选：常用 2 或 3；不限制就 0

# ============================================================
# 1) taxonomy: MAG_formal -> species_final
# ============================================================
tax <- readxl::read_xlsx(tax_path)

need_tax <- c("MAG_formal", "Species")
miss_tax <- setdiff(need_tax, colnames(tax))
if (length(miss_tax) > 0) stop("Taxonomy file missing columns: ", paste(miss_tax, collapse = ", "))

tax2 <- tax %>%
  transmute(
    MAG_formal = as.character(MAG_formal),
    species_final = as.character(Species)
  )

# ============================================================
# 2) LEfSe results
# ============================================================
lefse <- readxl::read_xlsx(lefse_path)

need_lefse <- c("MAG_formal", "Rumenotype", "LDA_score", "Pvalue")
miss_lefse <- setdiff(need_lefse, colnames(lefse))
if (length(miss_lefse) > 0) stop("LEfSe file missing columns: ", paste(miss_lefse, collapse = ", "))

lefse2 <- lefse %>%
  transmute(
    MAG_formal = as.character(MAG_formal),
    Rumenotype = as.character(Rumenotype),
    LDA_score  = suppressWarnings(as.numeric(LDA_score)),
    Pvalue = suppressWarnings(as.numeric(Pvalue))
    ) %>%
  filter(Rumenotype %in% c("Acetate_type", "Propionate_type")) %>%
  mutate(
    Rumenotype = factor(Rumenotype, levels = c("Acetate_type", "Propionate_type")),
    sig = LDA_score >= lda_cut
  )

# 只取显著富集的 MAG
lefse_sig <- lefse2 %>% filter(sig)

# ============================================================
# 3) Join by MAG_formal -> species
# ============================================================
dat <- lefse_sig %>%
  left_join(tax2, by = "MAG_formal") %>%
  mutate(
    species_final = ifelse(is.na(species_final) | species_final == "", MAG_formal, species_final),
    # 让 Propionate_type 朝下（负值），Acetate_type 朝上（正值）
    LDA_signed = ifelse(Rumenotype == "Propionate_type", -abs(LDA_score), abs(LDA_score))
  )

# 处理重复 species 标签：做成唯一标签（否则 x 轴可能重复导致显示/聚合问题）
dat <- dat %>%
  mutate(
    species_label = paste0(species_final, " (", MAG_formal, ")")
  )




# ============================================================
# 4) 排序规则（关键）
# - Acetate_type：LDA 大的在左侧（从大到小）
# - Propionate_type：LDA 大的在右侧
#   => 为了“越大越靠右”，Propionate 的因子水平要“从小到大”
# - 整体顺序：先排 Acetate，再排 Propionate
# ============================================================
levels_acet <- dat %>%
  filter(Rumenotype == "Acetate_type") %>%
  arrange(desc(LDA_score)) %>%
  pull(species_label)

levels_prop <- dat %>%
  filter(Rumenotype == "Propionate_type") %>%
  arrange(LDA_score) %>%             # 小 -> 大，确保大的在右侧
  pull(species_label)

dat <- dat %>%
  mutate(
    species_label = factor(species_label, levels = c(levels_acet, levels_prop))
  )

# ============================================================
# 5) Plot
# ============================================================
p <- ggplot(dat %>% subset(!str_detect(species_final, 'Others')), 
            aes(y = species_label, x = LDA_signed, fill = Rumenotype)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = type_cols) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    x = "LDA score",
    fill = "Rumenotype"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    axis.title.y = element_blank()
  )

print(p)

# 可选保存
ggsave("./cluster_meanRA_tax_ranks_bar/LEfSe_LDA_by_Species_barplot.pdf", p, width = 8, height = 6, dpi = 300)
write.csv(dat, "LEfSe_LDA_by_Species_for_plot.csv", row.names = FALSE)

dat <- left_join(dat, tax, by= 'MAG_formal')

