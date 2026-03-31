#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# =========================
# 0) 输入/输出（按需修改）
# =========================
in_xlsx <- "./FinalFigure/Supplemental file 4-AEE PEE differential MAGs.xlsx"
sheet_name <- 1   # 也可以写成 "Sheet1"
out_tsv  <- "./微生物/diffMAG_count_by_rumenotype_cluster.tsv"
out_pdf  <- "./微生物/figures/diffMAG_bubble_by_rumenotype_cluster.pdf"
out_png  <- "./微生物/figures/diffMAG_bubble_by_rumenotype_cluster.png"

# 你期望的 cluster 范围
clusters_keep <- 1:16

# =========================
# 1) 读入数据
# =========================
df <- readxl::read_xlsx(in_xlsx, sheet = sheet_name)

# 兼容列名（你表里一般会有这些）
# 必需列：Rumenotype, MAG.cluster, MAG_formal(或 MAG)
stopifnot(any(names(df) %in% c("Rumenotype")))
stopifnot(any(names(df) %in% c("MAG.cluster", "MAG_cluster", "cluster", "Cluster")))
stopifnot(any(names(df) %in% c("MAG_formal", "MAG", "mag", "MAG_ID", "MAG.id")))


df2 <- df %>%
  rename(
    MAG_cluster = `MAG.cluster`,
    MAG_id      = MAG_formal
  ) %>%
  mutate(
    Rumenotype  = as.character(Rumenotype),
    MAG_cluster = as.integer(MAG_cluster),
    MAG_id      = as.character(MAG_id)
  ) %>%
  filter(
    !is.na(Rumenotype),
    !is.na(MAG_cluster),
    !is.na(MAG_id),
    MAG_cluster %in% clusters_keep
  )

# =========================
# 2) 统计：每个 Rumenotype × cluster 的差异 MAG 数量（distinct MAG）
# =========================
cnt <- df2 %>%
  distinct(Rumenotype, MAG_cluster, MAG_id) %>%
  count(Rumenotype, MAG_cluster, name = "n_diffMAG") %>%
  tidyr::complete(
    Rumenotype,
    MAG_cluster = clusters_keep,
    fill = list(n_diffMAG = 0)
  )

# 导出表
readr::write_tsv(cnt, out_tsv)

# =========================
# 3) 气泡图（X=Rumenotype, Y=cluster 1-16, size=数量）
#    让 cluster=1 在最上面：y 设为 factor(levels=16:1)
# =========================
plot_df <- cnt %>%
  mutate(
    MAG_cluster_f = factor(MAG_cluster, levels = rev(clusters_keep)),
    Rumenotype_f  = factor(Rumenotype, levels = sort(unique(Rumenotype)))
  )

p <- ggplot(plot_df, aes(x = Rumenotype_f, y = MAG_cluster_f, size = n_diffMAG)) +
  geom_point(alpha = 0.8) +
  scale_size_area(max_size = 14) +
  labs(
    x = "Rumenotype",
    y = "MAG.cluster (1–16)",
    size = "Diff MAG count"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p

ggsave(out_pdf, p, width = 6.5, height = 6.5, useDingbats = FALSE)
ggsave(out_png, p, width = 6.5, height = 6.5, dpi = 300)

message("Done!\n- table: ", out_tsv, "\n- plot: ", out_pdf, " / ", out_png)
