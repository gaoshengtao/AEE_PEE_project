# ============================================================
# Project: Project LEfSe MAGs -> tSNE space
# Author: (you)
# ============================================================

# ---- packages ----
pkgs <- c("readr", "readxl", "dplyr", "ggplot2", "stringr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# ---- file paths ----
emb_path  <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/MAG_embeddings_clusters.tsv"
lefse_path <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/lefse.LDA_repMAGs_tpm.xlsx"

# ---- colors (as you provided) ----
type_cols <- c(
  "Acetate_type"    = "#2F5597",
  "Propionate_type" = "#C00000"
)

# ============================================================
# 1) Read embeddings (tsv): get TSNE1, TSNE2
# ============================================================
emb <- readr::read_tsv(emb_path, show_col_types = FALSE)

need_cols_emb <- c("MAG", "TSNE1", "TSNE2")
miss_emb <- setdiff(need_cols_emb, colnames(emb))
if (length(miss_emb) > 0) {
  stop("Embeddings table is missing columns: ", paste(miss_emb, collapse = ", "))
}

emb2 <- emb %>%
  transmute(
    MAG = as.character(MAG),
    MAG_formal = if ("MAG_formal" %in% colnames(emb)) as.character(MAG_formal) else NA_character_,
    TSNE1 = as.numeric(TSNE1),
    TSNE2 = as.numeric(TSNE2)
  )

# ============================================================
# 2) Read LEfSe results (xlsx)
# ============================================================
lefse <- readxl::read_xlsx(lefse_path)

need_cols_lefse <- c("MAG", "Rumenotype", "LDA_score", "Pvalue")
miss_lefse <- setdiff(need_cols_lefse, colnames(lefse))
if (length(miss_lefse) > 0) {
  stop("LEfSe table is missing columns: ", paste(miss_lefse, collapse = ", "))
}

alpha_p <- 0.05
lda_cut <- 0    # 如果你希望加 LDA 阈值，改成 2 或 3（常用 2）

lefse2 <- lefse %>%
  filter(Rumenotype %in% c("Acetate_type", "Propionate_type")) %>%
  transmute(
    MAG = as.character(MAG),
    Rumenotype = as.character(Rumenotype),
    LDA_score = as.numeric(LDA_score),
    pvalue = as.numeric(Pvalue)
  ) %>%
  mutate(
    Rumenotype = factor(Rumenotype, levels = c("Acetate_type", "Propionate_type")),
    sig = LDA_score >= lda_cut
  )



# ============================================================
# 3) Join: project LEfSe MAGs onto TSNE space
#    - First try join by MAG
#    - If too many missing, try join by MAG_formal
# ============================================================
# 先把 lefse 结果合并到所有 emb2（保持 emb2 的所有 MAG）
dat_all <- emb2 %>%
  left_join(lefse2, by = "MAG")

# 如果 join 缺失太多，尝试 MAG_formal
missing_rate_mag <- mean(is.na(dat_all$Rumenotype))
if (missing_rate_mag > 0.3 && any(!is.na(emb2$MAG_formal))) {
  message(sprintf("Join by MAG missing LEfSe rate = %.1f%%; trying MAG_formal...", 100 * missing_rate_mag))
  
  emb_alt <- emb2 %>%
    filter(!is.na(MAG_formal) & MAG_formal != "") %>%
    select(MAG_formal, MAG, TSNE1, TSNE2)
  
  dat_all2 <- emb_alt %>%
    left_join(lefse2, by = c("MAG_formal" = "MAG"))
  
  if (mean(is.na(dat_all2$Rumenotype)) < missing_rate_mag) {
    dat_all <- dat_all2 %>%
      transmute(
        MAG = MAG,
        TSNE1 = TSNE1,
        TSNE2 = TSNE2,
        Rumenotype = Rumenotype,
        LDA_score = LDA_score,
        Pvalue = Pvalue,
        sig = sig
      )
    message("Using MAG_formal join.")
  } else {
    message("MAG_formal join not better; keep join by MAG.")
  }
}

# 只保留有 tSNE 坐标的
dat_all <- dat_all %>% filter(!is.na(TSNE1), !is.na(TSNE2))

# 显著点数据
dat_sig <- dat_all %>% filter(sig == TRUE)

taxon <- openxlsx::read.xlsx("./master_repMAG_taxonomy_by_rank.xlsx")

dat_sig <- left_join(dat_sig, taxon, by = 'MAG_formal')

dat_sig$label <- ifelse(dat_sig$Species != 'Others', dat_sig$MAG_formal, NA) 

openxlsx::write.xlsx(dat_sig,'./AEE PEE differential MAGs.xlsx')


# ============================================================
# 4) Plot: LEfSe-enriched MAGs on TSNE1/TSNE2
#    - point size can reflect |LDA_score| (optional)
# ============================================================
p <- ggplot() +
  # 底图：所有 MAG（浅灰）
  geom_point(
    data = dat_all,
    aes(x = TSNE1, y = TSNE2),
    color = "grey80",
    size = 1.2,
    alpha = 0.6
  ) +
  # 叠加：显著 LEfSe MAG（彩色）
  geom_point(
    data = dat_sig,
    aes(x = TSNE1, y = TSNE2, color = Rumenotype),
    alpha = 0.95,
    size = 2
  ) +
  # ✅ 标签：只用 dat_sig 且 label 非 NA
  ggrepel::geom_text_repel(
    data = dat_sig %>% filter(!is.na(label)),
    aes(x = TSNE1, y = TSNE2, label = label, color = Rumenotype),
    size = 3,
    max.overlaps = Inf # 标签多的时候不让它“自动放弃”
  ) +
  scale_color_manual(values = type_cols) +
  labs(
    x = "t-SNE1",
    y = "t-SNE2",
    color = "Rumenotype"   # ✅ 这里用 color，不是 fill
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

print(p)


ggsave('./微生物/figures/lefse_MAGs_on_TSNE.pdf', p, width = 10, height = 10, dpi = 300)
# ============================================================
# 5) Save outputs (optional)
# ============================================================



# also save the merged table used for plotting
out_tsv <- "lefse_MAGs_on_TSNE.merged.tsv"
readr::write_tsv(dat_plot, out_tsv)

message("Done!\nSaved: ", out_png, "\nSaved: ", out_tsv)




