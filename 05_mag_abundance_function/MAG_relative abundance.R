#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
})

## =========================
## 0) Inputs (EDIT)
## =========================
ABUND_TSV   <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/repMAGs_coverm.tsv"
CLUSTER_TSV <- "./微生物/MAG_COG_cluster/Tables/MAG_embeddings_clusters.tsv"

REP_ANN_TSV <- "./微生物/drep_reps_phylophlan_annotation.tsv"
CDB_CSV     <- "./微生物/Cdb-dRep clusters.csv"
GTDB_BAC    <- "./微生物/gtdbtk.bac120.summary.tsv"
GTDB_AR     <- "./微生物/gtdbtk.ar53.summary.tsv"

OUTDIR <- "./cluster_meanRA_tax_ranks_bar"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

TOPN_TAXA <- 20
NOISE_KEYS <- c("0","Noise","noise","NOISE","NA","Unassigned","unassigned")

RANKS <- c("c","o","f","g")  # class/order/family/genus
RANK_LABEL <- c(c="Class", o="Order", f="Family", g="Genus")

## =========================
## 1) Helpers
## =========================
to_gid <- function(x){
  x <- basename(x)
  x <- sub("\\.(fa|fna|fasta|fas)$", "", x, ignore.case = TRUE)
  x
}

get_rank <- function(tax, rank = c("d","p","c","o","f","g","s")){
  rank <- match.arg(rank)
  if (is.na(tax) || tax == "" || grepl("Unclassified", tax, ignore.case = TRUE)) return(NA_character_)
  pat <- paste0(rank, "__([^;]+)")
  m <- stringr::str_match(tax, pat)
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

read_tsv_guess <- function(path){
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

# rank 名太怪/缺失时的兜底
clean_taxon <- function(x){
  x <- as.character(x)
  x <- ifelse(is.na(x) | x=="" | x=="NA", NA_character_, x)
  x
}

## =========================
## 2) Read cluster labels (repMAG -> cluster)
## =========================
clu <- read_tsv_guess(CLUSTER_TSV) %>% clean_names()

col_tip <- intersect(names(clu), c("tip","genome","mag","root_id","id"))[1]
col_cluster <- intersect(names(clu), c("cluster","hdbscan_cluster","functional_cluster","cluster_id"))[1]
if (is.na(col_tip) || is.na(col_cluster)) {
  stop("CLUSTER_TSV 必须包含 tip(或 genome/mag/root_id) 和 cluster(或 hdbscan_cluster) 两列。")
}

clu <- clu %>%
  transmute(
    tip = to_gid(.data[[col_tip]]),
    cluster_raw = as.character(.data[[col_cluster]])
  ) %>%
  mutate(cluster = ifelse(cluster_raw %in% NOISE_KEYS, "0", cluster_raw)) %>%
  distinct(tip, cluster)

message("[cluster] repMAG labeled = ", nrow(clu))

## =========================
## 3) Build taxonomy for each repMAG:
##    inherited (rep) > siblings majority (same dRep primary_cluster)
##    for ranks c/o/f/g
## =========================
# rep inherited taxonomy (already curated)
rep_ann <- read_tsv_guess(REP_ANN_TSV) %>%
  clean_names() %>%
  mutate(
    root_id = to_gid(root_id),
    gtdb_taxonomy_inherited = as.character(gtdb_taxonomy_inherited),
    d_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "d"))),
    p_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "p"))),
    c_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "c"))),
    o_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "o"))),
    f_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "f"))),
    g_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "g"))),
    s_inh = clean_taxon(map_chr(gtdb_taxonomy_inherited, ~get_rank(.x, "s"))),
  ) %>%
  select(root_id,d_inh,p_inh, c_inh, o_inh, f_inh, g_inh, s_inh)

# dRep mapping
cdb <- readr::read_csv(CDB_CSV, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    genome_id = to_gid(as.character(genome)),
    primary_cluster = as.integer(primary_cluster)
  ) %>%
  filter(!is.na(primary_cluster))

rep2primary <- cdb %>%
  select(genome_id, primary_cluster) %>%
  distinct() %>%
  rename(tip = genome_id)

# GTDB taxonomy for all genomes (bac + ar)
gtdb_bac <- read_tsv_guess(GTDB_BAC) %>%
  clean_names() %>%
  mutate(
    genome_id = to_gid(as.character(user_genome)),
    classification = as.character(classification),
    d_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "d"))),
    p_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "p"))),
    c_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "c"))),
    o_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "o"))),
    f_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "f"))),
    g_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "g"))),
    s_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "s")))
  ) %>%
  select(genome_id, d_gtdb, p_gtdb,c_gtdb, o_gtdb, f_gtdb, g_gtdb, s_gtdb)

gtdb_ar <- read_tsv_guess(GTDB_AR) %>%
  clean_names() %>%
  mutate(
    genome_id = to_gid(as.character(user_genome)),
    classification = as.character(classification),
    d_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "d"))),
    p_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "p"))),
    c_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "c"))),
    o_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "o"))),
    f_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "f"))),
    g_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "g"))),
    s_gtdb = clean_taxon(map_chr(classification, ~get_rank(.x, "s")))
  ) %>%
  select(genome_id, d_gtdb, p_gtdb,c_gtdb, o_gtdb, f_gtdb, g_gtdb, s_gtdb)

gtdb_all <- bind_rows(gtdb_bac, gtdb_ar) %>%
  distinct(genome_id, .keep_all = TRUE)

# siblings majority for each rank within primary_cluster
cdb_tax <- cdb %>%
  left_join(gtdb_all, by = "genome_id") %>%
  filter(!is.na(primary_cluster))

majority_rank <- function(df, rank_col, out_col){
  df %>%
    filter(!is.na(.data[[rank_col]]) & .data[[rank_col]] != "") %>%
    count(primary_cluster, !!sym(rank_col), sort = TRUE) %>%
    group_by(primary_cluster) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(primary_cluster, !!out_col := .data[[rank_col]])
}

maj_d <- majority_rank(cdb_tax, "d_gtdb", "d_sib")
maj_p <- majority_rank(cdb_tax, "p_gtdb", "p_sib")
maj_c <- majority_rank(cdb_tax, "c_gtdb", "c_sib")
maj_o <- majority_rank(cdb_tax, "o_gtdb", "o_sib")
maj_f <- majority_rank(cdb_tax, "f_gtdb", "f_sib")
maj_g <- majority_rank(cdb_tax, "g_gtdb", "g_sib")
maj_s <- majority_rank(cdb_tax, "s_gtdb", "s_sib")

# master taxonomy table for repMAGs
master_tax <- rep2primary %>%
  left_join(rep_ann, by = c("tip" = "root_id")) %>%
  left_join(maj_d, by = "primary_cluster") %>%
  left_join(maj_p, by = "primary_cluster") %>%
  left_join(maj_c, by = "primary_cluster") %>%
  left_join(maj_o, by = "primary_cluster") %>%
  left_join(maj_f, by = "primary_cluster") %>%
  left_join(maj_g, by = "primary_cluster") %>%
  left_join(maj_s, by = "primary_cluster") %>%
  left_join(clu, by = "tip") %>%
  mutate(
    domain_final  = ifelse(!is.na(d_inh), d_inh, ifelse(!is.na(d_sib), d_sib, "Others")),
    phylum_final  = ifelse(!is.na(p_inh), p_inh, ifelse(!is.na(p_sib), p_sib, "Others")),
    class_final  = ifelse(!is.na(c_inh), c_inh, ifelse(!is.na(c_sib), c_sib, "Others")),
    order_final  = ifelse(!is.na(o_inh), o_inh, ifelse(!is.na(o_sib), o_sib, "Others")),
    family_final = ifelse(!is.na(f_inh), f_inh, ifelse(!is.na(f_sib), f_sib, "Others")),
    genus_final  = ifelse(!is.na(g_inh), g_inh, ifelse(!is.na(g_sib), g_sib, "Others")),
    species_final  = ifelse(!is.na(s_inh), s_inh, ifelse(!is.na(s_sib), s_sib, "Others")),
  ) %>%
  filter(!is.na(cluster)) %>%
  select(tip, cluster, domain_final, phylum_final, class_final, order_final, family_final, genus_final, species_final) %>%
  distinct()

message("[master_tax] repMAG in plot = ", nrow(master_tax))

write_tsv(master_tax, file.path(OUTDIR, "master_repMAG_taxonomy_by_rank2.tsv"))

## =========================
## 4) Read CoverM RA% and compute mean RA% across samples for each repMAG
## =========================
abund <- read_tsv_guess(ABUND_TSV) %>% clean_names()

rep_col <- names(abund)[1]
abund <- abund %>% rename(tip = all_of(rep_col)) %>% mutate(tip = to_gid(tip))

ra_cols <- names(abund)[str_detect(names(abund), "relative_abundance")]
if (length(ra_cols) == 0) stop("未找到 Relative Abundance (%) 列（clean_names 后应包含 relative_abundance）。")

rep_mean_ra <- abund %>%
  select(tip, all_of(ra_cols)) %>%
  mutate(across(all_of(ra_cols), as.numeric)) %>%
  rowwise() %>%
  mutate(mean_ra_pct = mean(c_across(all_of(ra_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(tip, mean_ra_pct)

write_tsv(rep_mean_ra, file.path(OUTDIR, "repMAG_mean_relative_abundance.tsv"))

## =========================
## 5) Plot function (rank-wise)
## =========================
plot_rank <- function(rank_key){
  rank_col <- switch(rank_key,
                     c = "class_final",
                     o = "order_final",
                     f = "family_final",
                     g = "genus_final")
  rank_name <- RANK_LABEL[[rank_key]]
  
  # cluster x taxon: sum(repMAG mean RA%) within cluster
  df <- rep_mean_ra %>%
    inner_join(master_tax %>% select(tip, cluster, !!sym(rank_col)), by = "tip") %>%
    rename(taxon = !!sym(rank_col)) %>%
    group_by(cluster, taxon) %>%
    summarise(mean_ra_pct_sum = sum(mean_ra_pct, na.rm = TRUE), .groups = "drop")
  
  # pick top taxa globally (exclude Others first, then merge back)
  top_taxa <- df %>%
    filter(taxon != "Others") %>%
    group_by(taxon) %>%
    summarise(total_ra = sum(mean_ra_pct_sum, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_ra)) %>%
    slice_head(n = TOPN_TAXA) %>%
    pull(taxon)
  
  plot_df <- df %>%
    mutate(taxon_plot = ifelse(taxon %in% top_taxa, taxon, "Others")) %>%
    group_by(cluster, taxon_plot) %>%
    summarise(mean_ra_pct_sum = sum(mean_ra_pct_sum, na.rm = TRUE), .groups = "drop")
  
  # order taxa by global abundance (stack+legend order high->low)
  taxon_levels <- plot_df %>%
    group_by(taxon_plot) %>%
    summarise(total_ra = sum(mean_ra_pct_sum, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_ra)) %>%
    pull(taxon_plot)
  
  plot_df <- plot_df %>%
    mutate(taxon_plot = factor(taxon_plot, levels = taxon_levels))
  
  # order cluster: 1..N then 0 last, and reverse for coord_flip so 1 appears on top
  cluster_levels <- plot_df %>%
    distinct(cluster) %>%
    mutate(
      cluster_chr = as.character(cluster),
      is_noise = cluster_chr %in% c("0"),
      cluster_num = suppressWarnings(as.numeric(cluster_chr))
    ) %>%
    arrange(is_noise, cluster_num, cluster_chr) %>%
    pull(cluster_chr)
  
  plot_df <- plot_df %>%
    mutate(cluster = factor(as.character(cluster), levels = rev(cluster_levels)))
  
  # colors for top taxa (top20) + Others grey
  # 用固定的可重复调色板（20色），不够就循环
  base20 <- c(
    "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
    "#e6ab02","#a6761d","#1f78b4","#b2df8a","#fb9a99",
    "#cab2d6","#ffff99","#b15928","#8dd3c7","#bebada",
    "#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9"
  )
  
  # top taxa（不含 Others）按 abundance 顺序取颜色
  taxa_no_others <- setdiff(taxon_levels, "Others")
  cols <- rep(base20, length.out = length(taxa_no_others))
  names(cols) <- taxa_no_others
  cols <- c(cols, Others = "#BDBDBD")
  
  p <- ggplot(plot_df, aes(x = cluster, y = mean_ra_pct_sum, fill = taxon_plot)) +
    geom_col(width = 0.85) +
    coord_flip() +
    scale_fill_manual(
      values = cols,
      breaks = taxon_levels,   # legend 从高到低
      drop = FALSE
    ) +
    labs(
      x = "Cluster",
      y = "Mean relative abundance across samples (%)",
      fill = paste0("GTDB ", tolower(rank_name), " (top", TOPN_TAXA, ")"),
      title = paste0("Mean relative abundance of repMAGs in each cluster (stacked by ", tolower(rank_name), ")"),
      subtitle = "Per repMAG: mean(Relative Abundance %) across all samples; then summed within cluster"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  
  out_prefix <- file.path(OUTDIR, paste0("Fig_cluster_meanRA_by_", tolower(rank_name), "_top", TOPN_TAXA))
  ggsave(paste0(out_prefix, ".pdf"), p, width = 12, height = 6)
  ggsave(paste0(out_prefix, ".png"), p, width = 12, height = 6, dpi = 300)
  
  # output table
  write_tsv(plot_df %>% mutate(rank = rank_name),
            file.path(OUTDIR, paste0("cluster_meanRA_", tolower(rank_name), "_top", TOPN_TAXA, ".tsv")))
  
  message("[done] ", rank_name, " plot -> ", out_prefix, ".[pdf/png]")
}

## =========================
## 6) Run all ranks
## =========================
for (rk in RANKS) plot_rank(rk)

message("\nAll done. Outputs in: ", OUTDIR)
