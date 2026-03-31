#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

## =========================
## 0) Input / Output
## =========================
in_tsv <- "./微生物/MAG_COG_cluster/protein_similarity_with_taxonomy.full.tsv"
out_dir <- "protein_identity_density_by_class_pdfs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 只保留的 class（其余并入 Others）
keep_classes <- c(
  "Clostridia",
  "Bacteroidia",
  "Bacilli",
  "Methanobacteria",
  "Saccharimonadia",
  "Spirochaetia",
  "Alphaproteobacteria",
  "Coriobacteriia",
  "Actinomycetes",
  "Others"
)

## =========================
## 1) Read + clean
## =========================
dia_table <- readr::read_tsv(in_tsv, show_col_types = FALSE) %>%
  mutate(
    # 兼容 class_final / Class
    class_raw = case_when(
      "class_final" %in% names(.) ~ as.character(.data$Class),
      "Class" %in% names(.)       ~ as.character(.data$Class),
      TRUE                        ~ NA_character_
    ),
    MAG_formal = as.character(MAG_formal),
    pident     = as.numeric(pident)
  ) %>%
  filter(!is.na(class_raw), !is.na(MAG_formal), !is.na(pident)) %>%
  filter(pident >= 0, pident <= 100) %>%
  mutate(
    class_final = ifelse(class_raw %in% keep_classes, class_raw, "Others"),
    class_final = factor(class_final, levels = keep_classes)
  )

## 每个 MAG 至少 N 条蛋白（避免 density 失真）
min_n_per_mag <- 10
dia_table <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  filter(n() >= min_n_per_mag) %>%
  ungroup()

## =========================
## 2) 每个 MAG 的 mean pident（用于颜色映射）
## =========================
mag_color <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  summarise(
    pident_mean = mean(pident, na.rm = TRUE),
    .groups = "drop"
  )

## =========================
## 3) Density + Relative density
## =========================
dens_df <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  summarise(
    dens = list(density(pident, from = 0, to = 100, n = 200, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(dens, ~ .x$x),
    y = map(dens, ~ .x$y)
  ) %>%
  select(-dens) %>%
  unnest(c(x, y)) %>%
  group_by(class_final, MAG_formal) %>%
  mutate(rel_y = y / max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mag_color, by = c("class_final", "MAG_formal"))

## =========================
## 4) One class -> one PDF
##    颜色：浅蓝 → 浅红（按 pident_mean）
## =========================
for (cls in keep_classes) {
  
  dfc <- dens_df %>% filter(class_final == cls)
  if (nrow(dfc) == 0) {
    message("[WARN] No data for class: ", cls, " (skip)")
    next
  }
  
  p <- ggplot(
    dfc,
    aes(
      x = x,
      y = rel_y,
      group = MAG_formal,
      color = pident_mean
    )
  ) +
    geom_line(linewidth = 0.5, alpha = 0.9) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    scale_color_gradientn(
      colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
      limits = c(75, 100),
      oob = scales::squish,
      name = "Mean protein identity (%)"
    ) +
    labs(
      title = cls,
      x = "Protein identity (%)",
      y = "Relative density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  out_pdf <- file.path(
    out_dir,
    paste0("protein_identity_density_", gsub("[^A-Za-z0-9_]+", "_", cls), ".pdf")
  )
  
  ggsave(out_pdf, p, width = 3.6, height = 3.6, device = cairo_pdf)
  message("[OK] Saved: ", out_pdf)
}

message("[DONE] All class PDFs in: ", out_dir)#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
})

## =========================
## 0) Input / Output
## =========================
in_tsv <- "./微生物/MAG_COG_cluster/protein_similarity_with_taxonomy.full.tsv"
out_dir <- "protein_identity_density_by_class_pdfs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 只保留的 class（其余并入 Others）
keep_classes <- c(
  "Clostridia",
  "Bacteroidia",
  "Bacilli",
  "Methanobacteria",
  "Saccharimonadia",
  "Spirochaetia",
  "Alphaproteobacteria",
  "Coriobacteriia",
  "Actinomycetes",
  "Others"
)

## =========================
## 1) Read + clean
## =========================
dia_table <- readr::read_tsv(in_tsv, show_col_types = FALSE) %>%
  mutate(
    # 兼容 class_final / Class
    class_raw = case_when(
      "class_final" %in% names(.) ~ as.character(.data$class_final),
      "Class" %in% names(.)       ~ as.character(.data$Class),
      TRUE                        ~ NA_character_
    ),
    MAG_formal = as.character(MAG_formal),
    pident     = as.numeric(pident)
  ) %>%
  filter(!is.na(class_raw), !is.na(MAG_formal), !is.na(pident)) %>%
  filter(pident >= 0, pident <= 100) %>%
  mutate(
    class_final = ifelse(class_raw %in% keep_classes, class_raw, "Others"),
    class_final = factor(class_final, levels = keep_classes)
  )

## 每个 MAG 至少 N 条蛋白（避免 density 失真）
min_n_per_mag <- 10
dia_table <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  filter(n() >= min_n_per_mag) %>%
  ungroup()

## =========================
## 2) 每个 MAG 的 mean pident（用于颜色映射）
## =========================
mag_color <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  summarise(
    pident_mean = mean(pident, na.rm = TRUE),
    .groups = "drop"
  )

## =========================
## 3) Density + Relative density
## =========================
dens_df <- dia_table %>%
  group_by(class_final, MAG_formal) %>%
  summarise(
    dens = list(density(pident, from = 0, to = 100, n = 200, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(dens, ~ .x$x),
    y = map(dens, ~ .x$y)
  ) %>%
  select(-dens) %>%
  unnest(c(x, y)) %>%
  group_by(class_final, MAG_formal) %>%
  mutate(rel_y = y / max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(mag_color, by = c("class_final", "MAG_formal"))

## =========================
## 4) One class -> one PDF
##    颜色：浅蓝 → 浅红（按 pident_mean）
## =========================
for (cls in keep_classes) {
  
  dfc <- dens_df %>% filter(class_final == cls)
  if (nrow(dfc) == 0) {
    message("[WARN] No data for class: ", cls, " (skip)")
    next
  }
  
  p <- ggplot(
    dfc,
    aes(
      x = x,
      y = rel_y,
      group = MAG_formal,
      color = pident_mean
    )
  ) +
    geom_line(linewidth = 0.5, alpha = 0.9) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    scale_color_gradientn(
      colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
      limits = c(75, 100),
      oob = scales::squish,
      name = "Mean protein identity (%)"
    ) +
    labs(
      title = cls,
      x = "Protein identity (%)",
      y = "Relative density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  out_pdf <- file.path(
    out_dir,
    paste0("protein_identity_density_", gsub("[^A-Za-z0-9_]+", "_", cls), ".pdf")
  )
  
  ggsave(out_pdf, p, width = 3, height = 3, device = cairo_pdf)
  message("[OK] Saved: ", out_pdf)
}

message("[DONE] All class PDFs in: ", out_dir)