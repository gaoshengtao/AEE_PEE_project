suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

strip_ext <- function(x){
  sub("\\.(fa|fna|fasta|fas|fsa)$", "", x, ignore.case = TRUE)
}

# -----------------------
# 1) Cluster vs rest: KO enrichment
# -----------------------
cluster_ko_enrichment <- function(
    ko_mat_tsv,
    embed_tsv,
    cluster_col = "hdbscan_cluster",
    out_tsv = "Cluster_KO_enrichment.tsv",
    use_log1p = TRUE
){
  ko <- fread(ko_mat_tsv, data.table = FALSE)
  colnames(ko)[1] <- "MAG"
  ko$MAG <- strip_ext(as.character(ko$MAG))
  
  emb <- fread(embed_tsv, data.table = FALSE)
  emb$MAG <- strip_ext(as.character(emb$MAG))
  stopifnot(cluster_col %in% colnames(emb))
  
  dat <- inner_join(emb[, c("MAG", cluster_col)], ko, by="MAG")
  dat[[cluster_col]] <- as.character(dat[[cluster_col]])
  # 可选：把 "0" / "Noise" / NA 改名
  dat[[cluster_col]][dat[[cluster_col]] %in% c("0","Noise")] <- "NA"
  dat[[cluster_col]][is.na(dat[[cluster_col]])] <- "NA"
  
  feat_cols <- setdiff(colnames(dat), c("MAG", cluster_col))
  X <- as.matrix(dat[, feat_cols])
  storage.mode(X) <- "numeric"
  X[!is.finite(X)] <- 0
  
  if (use_log1p) {
    Xv <- log2(1 + X)  # 用 log2(1+copy) 做“表达量”
  } else {
    Xv <- X
  }
  
  cl <- dat[[cluster_col]]
  clusters <- sort(unique(cl))
  
  res_all <- list()
  
  for (cc in clusters) {
    in_idx  <- cl == cc
    out_idx <- cl != cc
    
    Xin  <- Xv[in_idx, , drop=FALSE]
    Xout <- Xv[out_idx, , drop=FALSE]
    
    mean_in  <- colMeans(Xin)
    mean_out <- colMeans(Xout)
    prev_in  <- colMeans(X[in_idx, , drop=FALSE] > 0)   # prevalence 用原始copy>0
    prev_out <- colMeans(X[out_idx, , drop=FALSE] > 0)
    
    # Wilcoxon (每个KO)
    pvals <- sapply(seq_len(ncol(Xv)), function(j){
      suppressWarnings(wilcox.test(Xv[in_idx, j], Xv[out_idx, j])$p.value)
    })
    
    df <- data.frame(
      cluster = cc,
      KO = colnames(Xv),
      mean_in = mean_in,
      prev_in = prev_in,
      mean_out = mean_out,
      prev_out = prev_out,
      log2FC = mean_in - mean_out, # 因为是log2(1+copy)，直接差值就是log2FC
      pvalue = pvals,
      stringsAsFactors = FALSE
    )
    df$fdr <- p.adjust(df$pvalue, method="BH")
    
    res_all[[cc]] <- df
  }
  
  res <- bind_rows(res_all)
  fwrite(res, out_tsv, sep="\t")
  res
}

# -----------------------
# 2) Module enrichment from significant KO
# -----------------------
cluster_module_enrichment_fromKO <- function(
    ko_enrich_df,
    module_ko_map,
    module_name_file = "微生物/MAG_COG_cluster/Tables/module name.txt",
    log2fc_cut = 1,
    fdr_cut = 0.05,
    min_k_overlap = 2,
    min_module_size = 3,
    fallback_top_ko = 50,
    topN_modules_per_cluster = 10,
    out_tsv = NULL
){
  suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
  })
  
  # ---------- 0) detect columns (robust) ----------
  cn <- colnames(ko_enrich_df)
  
  pick1 <- function(cands){
    hit <- intersect(cands, cn)
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  col_cluster <- pick1(c("cluster","Cluster","cluster_id","ClusterID"))
  col_ko      <- pick1(c("KO","ko","kegg_ko","KEGG_ko","term","id","feature"))
  col_log2fc  <- pick1(c("log2FC","avg_log2FC","avg_log2fc","logFC","log2fc"))
  col_fdr     <- pick1(c("fdr","FDR","padj","qvalue","q_value"))
  
  message("[detect columns] cluster=", col_cluster,
          "; KO=", col_ko,
          "; log2FC=", col_log2fc,
          "; fdr=", col_fdr)
  
  if (any(is.na(c(col_cluster, col_ko, col_log2fc, col_fdr)))) {
    stop("ko_enrich_df 缺少必要列。请先看：colnames(ko_enrich_df)")
  }
  
  # ---------- 1) KO enrich table -> DT ----------
  df <- data.table(
    cluster = as.character(ko_enrich_df[[col_cluster]]),
    KO      = as.character(ko_enrich_df[[col_ko]]),
    log2FC  = suppressWarnings(as.numeric(ko_enrich_df[[col_log2fc]])),
    fdr     = suppressWarnings(as.numeric(ko_enrich_df[[col_fdr]]))
  )
  df[, KO := str_replace_all(KO, "^ko:", "")]
  df[, KO := str_trim(KO)]
  df <- df[!is.na(KO) & KO != "" & !is.na(log2FC) & !is.na(fdr)]
  
  # ---------- 2) module<->KO map ----------
  m <- fread(module_ko_map, data.table = TRUE)
  if (!all(c("Module","KO") %in% names(m))) {
    stop("module_ko_map.txt 必须包含列名：Module 和 KO")
  }
  
  m[, Module := as.character(Module)]
  m[, KO := as.character(KO)]
  m[, KO := str_replace_all(KO, "^ko:", "")]
  # 允许第二列里是 “K00001,K00002;K00003” 这种
  m[, KO := str_replace_all(KO, "[;|\\s]+", ",")]
  m[, KO := str_replace_all(KO, ",+", ",")]
  m[, KO := str_trim(KO)]
  
  m2 <- m[, .(KO = unlist(strsplit(KO, ",", fixed = TRUE))), by = Module]
  m2[, KO := str_trim(KO)]
  m2 <- unique(m2[!is.na(KO) & KO != ""])
  
  module_size <- m2[, .(K = .N), by = Module]
  module_size <- module_size[K >= min_module_size]
  
  m2 <- merge(m2, module_size, by = "Module", all = FALSE)
  
  bg_kos <- unique(m2$KO)
  N_bg <- length(bg_kos)
  
  clusters <- sort(unique(df$cluster))
  
  # ---------- 3) optional module names ----------
  mn <- NULL
  if (!is.null(module_name_file) && file.exists(module_name_file)) {
    mn <- fread(module_name_file, data.table = TRUE)
    # 兼容大小写
    if ("Name" %in% names(mn) && !"name" %in% names(mn)) setnames(mn, "Name", "name")
    if ("module" %in% names(mn) && !"Module" %in% names(mn)) setnames(mn, "module", "Module")
    if (!all(c("Module","name") %in% names(mn))) {
      mn <- NULL
    } else {
      mn[, Module := as.character(Module)]
      mn[, name := as.character(name)]
      mn <- unique(mn[, .(Module, name)])
    }
  }
  
  # ---------- 4) per-cluster enrichment ----------
  out_all <- vector("list", length(clusters))
  names(out_all) <- clusters
  
  for (cl in clusters) {
    dfc <- df[cluster == cl]
    
    # 4.1 pick KO list: sig first, else fallback top log2FC
    sig <- dfc[log2FC >= log2fc_cut & fdr <= fdr_cut][order(-log2FC)]
    kos_use <- sig$KO
    kos_use <- kos_use[kos_use %in% bg_kos]
    mode <- "sig"
    
    if (length(unique(kos_use)) < min_k_overlap) {
      fb <- dfc[order(-log2FC)][1:min(fallback_top_ko, .N)]
      kos_use <- fb$KO
      kos_use <- kos_use[kos_use %in% bg_kos]
      kos_use <- unique(kos_use)
      mode <- "fallback"
    } else {
      kos_use <- unique(kos_use)
    }
    
    n_sig <- length(kos_use)
    
    if (n_sig == 0) {
      out_all[[cl]] <- data.table(
        cluster = cl, Module = NA_character_, name = NA_character_,
        Module_id_name = NA_character_, K = NA_integer_,
        n_sig = 0L, k = 0L, enrichment_ratio = NA_real_,
        pvalue = NA_real_, fdr = NA_real_, mode = mode
      )
      next
    }
    
    # 4.2 overlap count per module
    ov <- m2[KO %in% kos_use, .(k = .N), by = Module]
    ov <- merge(ov, module_size, by = "Module", all.x = TRUE, all.y = FALSE)
    
    # 最少重叠数：优先 min_k_overlap，否则放宽到 >=1 以保证每个cluster都有module可输出
    ov2 <- ov[k >= min_k_overlap]
    min_k_used <- min_k_overlap
    if (nrow(ov2) == 0L) {
      ov2 <- ov[k >= 1]
      min_k_used <- 1
    }
    
    if (nrow(ov2) == 0L) {
      out_all[[cl]] <- data.table(
        cluster = cl, Module = NA_character_, name = NA_character_,
        Module_id_name = NA_character_, K = NA_integer_,
        n_sig = n_sig, k = 0L, enrichment_ratio = NA_real_,
        pvalue = NA_real_, fdr = NA_real_,
        mode = paste0(mode, ";min_k_used=", min_k_used)
      )
      next
    }
    
    # 4.3 hypergeometric p + enrichment ratio
    N <- N_bg
    n <- n_sig
    
    ov2[, pvalue := phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)]
    ov2[, enrichment_ratio := (k / n) / (K / N)]
    ov2[, cluster := cl]
    ov2[, n_sig := n]
    ov2[, mode := paste0(mode, ";min_k_used=", min_k_used)]
    
    setorder(ov2, pvalue, -enrichment_ratio)
    ov2[, fdr := p.adjust(pvalue, method = "BH")]
    
    # 4.4 add module name
    if (!is.null(mn)) {
      ov2 <- merge(ov2, mn, by = "Module", all.x = TRUE)
    } else {
      ov2[, name := NA_character_]
    }
    ov2[is.na(name) | name == "", name := Module]
    ov2[, Module_id_name := paste0(Module, ": ", name)]
    
    # 4.5 guarantee output: topN per cluster even if not significant
    keepN <- min(topN_modules_per_cluster, nrow(ov2))
    out_all[[cl]] <- ov2[1:keepN, .(
      cluster, Module, name, Module_id_name, K, n_sig, k,
      enrichment_ratio, pvalue, fdr, mode
    )]
  }
  
  res <- rbindlist(out_all, use.names = TRUE, fill = TRUE)
  
  if (!is.null(out_tsv)) {
    dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
    fwrite(res, out_tsv, sep = "\t", quote = FALSE, na = "NA")
  }
  
  return(res)
}





ko_enrich <- cluster_ko_enrichment(
  ko_mat_tsv = "微生物/MAG_COG_cluster/Tables/MAGxFeature_copy_corrected.tsv",
  embed_tsv  = "微生物/MAG_COG_cluster/Tables/MAG_embeddings_clusters.tsv",
  cluster_col = "hdbscan_cluster",
  out_tsv = "微生物/MAG_COG_cluster/Tables/Cluster_KO_enrichment.tsv"
)



mod_enrich <- cluster_module_enrichment_fromKO(
  ko_enrich_df = ko_enrich,
  module_ko_map = "微生物/MAG_COG_cluster/Tables/module_ko_map.txt",
  module_name_file = "微生物/MAG_COG_cluster/Tables/module name.txt",
  log2fc_cut = 1,
  fdr_cut = 0.05,
  min_k_overlap = 2,
  min_module_size = 3,
  out_tsv = "微生物/MAG_COG_cluster/Tables/Cluster_Module_enrichment_fromKO.tsv"
)



