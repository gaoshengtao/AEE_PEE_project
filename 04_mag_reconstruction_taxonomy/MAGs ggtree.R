#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(readxl)
  library(ape)
  library(janitor)
  
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
})

## =========================
## 0) Inputs
## =========================
TREE_FILE <- "./微生物/tree_input_genomes_resolved.tre"
TAX_XLSX  <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 3-master_repMAG_taxonomy_by_rank.xlsx"
CDB_CSV   <- "./微生物/Cdb-dRep clusters.csv"
QUAL_CSV  <- "./微生物/genomeInformation.csv"
QUAL_RERUN_TSV <- "./微生物/MAG quality_report_rerun.tsv"
ABUND_TSV <- "./微生物/repMAGs_coverm.tsv"   # optional

OUT_PDF <- "./微生物/figures/ggtree_tracks.pdf"

## =========================
## 1) Helpers
## =========================
clean_tip_label <- function(x){
  x %>%
    as.character() %>%
    str_replace_all("^['\"]|['\"]$", "")
}

make_tip_key <- function(x){
  x %>%
    clean_tip_label() %>%
    str_replace("_dup\\d+$", "")
}

to_gid <- function(x){
  x <- basename(x)
  sub("\\.(fa|fna|fasta|fas)$", "", x, ignore.case = TRUE)
}

read_tsv_guess <- function(path){
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

first_non_na <- function(x){
  x2 <- x[!is.na(x) & x != ""]
  if (length(x2) == 0) NA else x2[1]
}

## =========================
## 2) Read tree
## =========================
tr <- read.tree(TREE_FILE)

tips_df <- tibble(
  tip_raw = tr$tip.label,
  tip     = clean_tip_label(tr$tip.label),   # ✅ 用于绘图对齐 label
  tip_key = make_tip_key(tr$tip.label)       # ✅ 用于 join 匹配
)

message("[tree] tips: ", nrow(tips_df))
message("[tree] unique tip: ", n_distinct(tips_df$tip))
message("[tree] unique tip_key: ", n_distinct(tips_df$tip_key))

## =========================
## 3) Taxonomy (Excel) -> keep only Class
## =========================
tax <- readxl::read_xlsx(TAX_XLSX) %>%
  clean_names() %>%
  mutate(
    tip_key = make_tip_key(as.character(original_magid)),
    class   = as.character(class)
  ) %>%
  group_by(tip_key) %>%
  summarise(class = first_non_na(class), .groups = "drop")

## =========================
## 4) dRep clusters -> cluster_size; rep2cluster unique by tip_key
## =========================
cdb <- readr::read_csv(CDB_CSV, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    genome_id = to_gid(as.character(genome)),
    tip_key   = make_tip_key(genome_id),
    primary_cluster = as.integer(primary_cluster)
  ) %>%
  filter(!is.na(primary_cluster))

cluster_size <- cdb %>%
  count(primary_cluster, name = "cluster_size")

# 关键：同一 tip_key 可能对应多个 primary_cluster（因为 tip_key 去 dup 了）
# 规则：取 cluster_size 最大的那一个（更合理）
rep2cluster <- cdb %>%
  distinct(tip_key, primary_cluster) %>%
  left_join(cluster_size, by = "primary_cluster") %>%
  arrange(tip_key, desc(cluster_size), primary_cluster) %>%
  group_by(tip_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(tip_key, primary_cluster)

## =========================
## 5) Quality (base + rerun fill) -> unique by tip_key
## =========================
qual <- readr::read_csv(QUAL_CSV, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    tip_key = make_tip_key(to_gid(as.character(mag))),
    completeness  = suppressWarnings(as.numeric(completeness)),
    contamination = suppressWarnings(as.numeric(contamination)),
    gc_content    = suppressWarnings(as.numeric(gc_content)),
    contig_n50    = suppressWarnings(as.numeric(contig_n50)),
    genome_size_bp= suppressWarnings(as.numeric(genome_size))
  ) %>%
  mutate(
    log10_n50 = ifelse(!is.na(contig_n50) & contig_n50 > 0, log10(contig_n50), NA_real_),
    genome_size_mb = ifelse(!is.na(genome_size_bp), genome_size_bp / 1e6, NA_real_)
  ) %>%
  select(tip_key, completeness, contamination, gc_content, contig_n50, log10_n50, genome_size_mb) %>%
  group_by(tip_key) %>%
  summarise(across(where(is.numeric), ~ if (all(is.na(.x))) NA_real_ else .x[which.max(!is.na(.x))]), .groups = "drop")

qual_rerun <- NULL
if (file.exists(QUAL_RERUN_TSV) && file.info(QUAL_RERUN_TSV)$size > 0) {
  qual_rerun <- readr::read_tsv(QUAL_RERUN_TSV, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      tip_key = make_tip_key(to_gid(as.character(name))),
      gc_content_rerun = suppressWarnings(as.numeric(gc_content)),
      contig_n50_rerun = suppressWarnings(as.numeric(contig_n50))
    ) %>%
    mutate(
      log10_n50_rerun = ifelse(!is.na(contig_n50_rerun) & contig_n50_rerun > 0,
                               log10(contig_n50_rerun), NA_real_)
    ) %>%
    select(tip_key, gc_content_rerun, contig_n50_rerun, log10_n50_rerun) %>%
    group_by(tip_key) %>%
    summarise(across(where(is.numeric),~ if (all(is.na(.x))) NA_real_ else .x[which.max(!is.na(.x))]),
              .groups="drop")
}
## =========================
## 6) Abundance stats (optional)
## =========================
has_abund <- file.exists(ABUND_TSV) && file.info(ABUND_TSV)$size > 0
ABUND_METRIC <- "relative_abundance_percent"  # or "tpm"
abund_stats <- NULL

if (has_abund) {
  abund_raw <- read_tsv_guess(ABUND_TSV) %>% clean_names()
  rep_col <- names(abund_raw)[1]
  abund_raw <- abund_raw %>% rename(tip = all_of(rep_col)) %>%
    mutate(tip = clean_tip_label(tip))
  
  metric_pattern <- paste0("_", ABUND_METRIC, "$")
  metric_cols <- names(abund_raw)[str_detect(names(abund_raw), metric_pattern)]
  if (length(metric_cols) == 0) stop("No columns matched metric: ", ABUND_METRIC)
  
  abund_mat <- abund_raw %>%
    select(tip, all_of(metric_cols)) %>%
    mutate(across(all_of(metric_cols), ~suppressWarnings(as.numeric(.x)))) %>%
    column_to_rownames("tip") %>%
    as.matrix()
  
  abund_stats <- tibble(
    tip = rownames(abund_mat),
    mean_abundance = as.numeric(rowMeans(abund_mat, na.rm = TRUE)),
    prevalence = as.numeric(rowMeans(abund_mat > 0, na.rm = TRUE))
  )
}

## =========================
## 7) Master (aligned to tree tips; no duplicated tips)
## =========================
class_colors <- c(
  "Clostridia"        = "#D73027",
  "Bacteroidia"       = "#4575B4",
  "Bacilli"           = "#1A9850",
  "Methanobacteria"   = "#7B3294",
  "Saccharimonadia"    = "#66C2A5",
  "Spirochaetia"      = "#F46D43",
  "Alphaproteobacteria"     = "#E78AC3",
  "Coriobacteriia"    = "#FDAE61",
  "Actinomycetes"     = "#A6611A",
  "Others"            = "#BDBDBD"
)



master <- tips_df %>%
  transmute(tip, tip_key) %>%
  left_join(tax, by = "tip_key") %>%
  left_join(rep2cluster, by = "tip_key") %>%
  left_join(cluster_size, by = "primary_cluster") %>%
  left_join(qual, by = "tip_key")

if (!is.null(qual_rerun)) {
  master <- master %>%
    left_join(qual_rerun, by = "tip_key") %>%
    mutate(
      gc_content = coalesce(gc_content, gc_content_rerun),
      contig_n50 = coalesce(contig_n50, contig_n50_rerun),
      log10_n50  = coalesce(log10_n50,  log10_n50_rerun)
    ) %>%
    select(-gc_content_rerun, -contig_n50_rerun, -log10_n50_rerun)
}

if (!is.null(abund_stats)) {
  master <- master %>% left_join(abund_stats, by = "tip")
}

master <- master %>%
  mutate(
    class_final = ifelse(!is.na(class) & class != "" & class %in% names(class_colors),
                         class, "Others"),
    class_color = unname(class_colors[class_final]),
    
    # completeness bins
    comp_bin = case_when(
      is.na(completeness) ~ "NA",
      completeness >= 90  ~ ">=90",
      completeness >= 50  ~ "50-90",
      TRUE ~ "<50"
    ),
    
    # cluster_size 4-bin
    cs_bin = case_when(
      is.na(cluster_size) ~ NA_character_,
      cluster_size >= 300 ~ "300",
      cluster_size >= 100 ~ "100",
      cluster_size >= 10  ~ "10",
      cluster_size >= 1   ~ "1",
      TRUE ~ NA_character_
    )
    )

master <- master %>%
  mutate(
    mean_abundance_trans = ifelse(mean_abundance > 0.5, 0.5, mean_abundance)
  )


# ✅ sanity: 不允许 master 出现重复 tip
if (nrow(master) != n_distinct(master$tip)) {
  dup <- master %>% count(tip) %>% filter(n > 1)
  print(dup)
  stop("master has duplicated tips -> join explosion. Fix upstream uniqueness.")
}
message("[master] rows: ", nrow(master))

## =========================
## 8) Plot: base tree
## =========================
p <- ggtree(tr, layout = "fan", size = 0.15, open.angle = 15)

## ---- (A) Tip-end circles: cluster_size (size bins) + fill=Class ----
tip_dat <- p$data %>%
  filter(isTip) %>%
  mutate(tip_key = make_tip_key(label)) %>%    # ✅ tree侧统一key
  left_join(
    master %>% transmute(tip_key, cs_bin, class_final),  # ✅ master侧用key
    by = "tip_key"
  )


# 你 tippoint 用的 scale_size_manual 的合法水平
# valid_cs <- c("1","10","100","300")
# valid_class <- names(class_colors)
# 
# bad_tip_dat <- tip_dat %>%
#   mutate(
#     bad_size  = is.na(cs_bin) | !(cs_bin %in% valid_cs),
#     bad_class = is.na(class_final) | !(class_final %in% valid_class)
#   ) %>%
#   filter(bad_size | bad_class) %>%
#   select(label, cs_bin, class_final, bad_size, bad_class)
# 
# bad_tip_dat
# nrow(bad_tip_dat)


p <- p +
  geom_tippoint(
    data = tip_dat,
    aes(size = cs_bin, fill = class_final),
    alpha = 0.95,
    shape = 21,
    color = "white",
    stroke = 0.15
  ) +
  scale_size_manual(
    values = c("1"=1.6, "10"=2.4, "100"=3.6, "300"=4.2),
    breaks = c("1","10","100","300"),
    name = "Cluster size"
  ) +
  scale_fill_manual(values = class_colors, name = "Class") +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color()

p

## ---- (B) Track 1: completeness bar (fill=bins) ----
p <- p +
  geom_fruit(
    data = master %>% rename(label = tip),
    geom = geom_col,
    mapping = aes(y = label, x = completeness, fill = comp_bin),
    orientation = "y",
    width = 0.55,
    offset = 0.02,
    pwidth = 0.08
  ) +
  scale_fill_manual(
    values = c("50-90"="#A6CEE3", ">=90"="#FEE08B", "<50"="#DDDDDD", "NA"="transparent"),
    breaks = c("50-90", ">=90"),
    name = "Completeness"
  ) +
  ggnewscale::new_scale_fill()

p

## ---- (D) Tracks 3-5: contamination / gc_content / log10_n50 as gradient tiles (same blue palette) ----

contam_low  <- "#F7F9FB"
contam_high <- "#90B4E0"
gc_low  <- "#F8F6FA"
gc_high <- "#6E5A7E"
n50_low  <- "#F6FAF7"
n50_high <- "#4F7C6B"

# contamination
p <- p +
  geom_fruit(
    data = master %>% rename(label = tip),
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = contamination),
    width = 0.02,
    offset = 0.015,
    pwidth = 0.015
  ) +
  scale_fill_gradient(low = contam_low, high = contam_high, na.value = "transparent", name = "Contam.") +
  ggnewscale::new_scale_fill()


# gc_content
p <- p +
  geom_fruit(
    data = master %>% rename(label = tip),
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = gc_content),
    width = 0.02,
    offset = 0.015,
    pwidth = 0.015
  ) +
  scale_fill_gradient(low = gc_low, high = gc_high, na.value = "transparent", name = "GC") +
  ggnewscale::new_scale_fill()

# log10_n50
p <- p +
  geom_fruit(
    data = master %>% rename(label = tip),
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = log10_n50),
    width = 0.02,
    offset = 0.015,
    pwidth = 0.015
  ) +
  scale_fill_gradient(low = n50_low, high = n50_high, na.value = "transparent", name = "log10N50") +
  ggnewscale::new_scale_fill()

p


## ---- (E) Tracks 6-7: mean_abundance / prevalence as bar (another color family; no log) ----
  green_bar <- "#2CA25F"
  cyan_bar  <- "#41B6C4"
  
  # mean abundance bar
  p <- p +
    geom_fruit(
      data = master %>% rename(label = tip),
      geom = geom_col,
      mapping = aes(y = label, x = mean_abundance_trans),
      orientation = "y",
      width = 0.55,
      offset = 0.02,
      pwidth = 0.05,
      fill = green_bar
    )
  
  # prevalence bar
  p <- p +
    geom_fruit(
      data = master %>% rename(label = tip),
      geom = geom_col,
      mapping = aes(y = label, x = prevalence),
      orientation = "y",
      width = 0.55,
      offset = 0.02,
      fill = cyan_bar,
      pwidth = 0.05
    )



## ---- (C) Track 2: genome_size(Mb) bar (fill=Class color) ----
p <- p +
  geom_fruit(
    data = master %>% rename(label = tip),
    geom = geom_col,
    mapping = aes(y = label, x = genome_size_mb, fill = class_final),
    orientation = "y",
    width = 0.55,
    offset = 0.02,
    pwidth = 0.2
  ) +
  scale_fill_manual(values = class_colors, name = "Class") +
  ggnewscale::new_scale_fill()

p
## =========================
## 9) Export
## =========================
ggsave(OUT_PDF, p, width = 12, height = 12, units = "in", dpi = 300)
message("[out] ", OUT_PDF)
