# =========================
# LEfSe input from CoverM + rumenotype assignments
# Author: (you)
# =========================

# ---- packages ----
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(readr)
})

# ---- input paths (Windows) ----
coverm_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/repMAGs_coverm.tsv"
assign_xlsx <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/VFA_rumenotype_assignments.xlsx"

# ---- output path ----
out_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/lefse_input_repMAGs_RelAb.tsv"

# =========================
# 1) Read group assignments
# =========================
meta <- readxl::read_excel(assign_xlsx, sheet = 1) %>%
  as.data.frame()

# keep needed columns
# meta columns: SampleID k ClusterID Rumenotype TVFA AP ABP BCVFA
meta <- meta %>%
  dplyr::select(SampleID, Rumenotype) %>%
  mutate(
    SampleID = as.character(SampleID),
    # extract digits for mapping: D0619 -> 0619
    digits = str_extract(SampleID, "\\d+"),
    # keep only target classes
    Rumenotype = as.character(Rumenotype)
  ) %>%
  filter(Rumenotype %in% c("Acetate-type", "Propionate-type"))


if (nrow(meta) == 0) {
  stop("在 Sheet1 中没有找到 Rumenotype 为“乙酸型/丙酸型”的样本，请检查拼写是否一致。")
}

# =========================
# 2) Read CoverM table
# =========================
# coverm header example:
# Genome  H0619.repMAGs.sorted Mean  ...  H0619.repMAGs.sorted Relative Abundance (%) ...
coverm <- readr::read_tsv(coverm_tsv, show_col_types = FALSE)

if (!"Genome" %in% colnames(coverm)) {
  stop("repMAGs_coverm.tsv 中未找到 Genome 列，请检查文件格式/表头。")
}

# =========================
# 3) Choose abundance column for LEfSe
# =========================
# Choose Relative Abundance (%) columns (recommended for LEfSe)
abund_suffix <- "TPM" #"Relative Abundance (%)"

abund_cols <- colnames(coverm)[str_detect(colnames(coverm), fixed(abund_suffix, ignore_case = FALSE))]

if (length(abund_cols) == 0) {
  stop("在 repMAGs_coverm.tsv 中未找到 'Relative Abundance (%)' 列。请确认列名是否完全一致。")
}

# Parse sample code from column name, e.g.:
# "H0619.repMAGs.sorted Relative Abundance (%)" -> "H0619"
sample_codes <- abund_cols %>%
  str_replace("\\.repMAGs\\.sorted.*$", "")  # remove from ".repMAGs.sorted ..." to end

# Extract digits from sample code: H0619 -> 0619
sample_digits <- str_extract(sample_codes, "\\d+")

# Build mapping table for coverm columns
map_coverm <- data.frame(
  coverm_col = abund_cols,
  coverm_code = sample_codes,
  digits = sample_digits,
  stringsAsFactors = FALSE
)

# =========================
# 4) Match samples by digits
# =========================
# meta: SampleID like D0619, digits=0619
# coverm: digits from H0619, digits=0619
matched <- meta %>%
  inner_join(map_coverm, by = "digits")

if (nrow(matched) == 0) {
  stop("没有匹配到任何样本：请确认 repMAGs_coverm.tsv 列名里的样本号(如 H0619)与 meta 的 SampleID (如 D0619)数字部分一致。")
}

# Optional: warn about unmatched
unmatched_meta <- meta %>% anti_join(matched %>% distinct(SampleID, digits), by = c("SampleID", "digits"))
unmatched_cov  <- map_coverm %>% anti_join(matched %>% distinct(coverm_col, digits), by = c("coverm_col", "digits"))

if (nrow(unmatched_meta) > 0) {
  message("⚠️ 有这些 meta 样本未在 coverm 中匹配到，将被忽略：\n  ",
          paste(unmatched_meta$SampleID, collapse = ", "))
}
if (nrow(unmatched_cov) > 0) {
  message("⚠️ 有这些 coverm 样本列未在 meta 中匹配到（或不属于乙酸型/丙酸型），将被忽略：\n  ",
          paste(unique(unmatched_cov$coverm_code), collapse = ", "))
}

# Keep order (you can change ordering if you want)
# Here: keep as meta order, within that by appearance
matched <- matched %>%
  mutate(meta_order = match(SampleID, meta$SampleID)) %>%
  arrange(meta_order)

# Final sample list for output
final_sample_ids <- matched$SampleID
final_class      <- matched$Rumenotype
final_coverm_cols <- matched$coverm_col

# =========================
# 5) Build LEfSe format
# =========================
# Feature table: rows=Genome(MAG), cols=samples
feat <- coverm %>%
  dplyr::select(Genome, all_of(final_coverm_cols))

# In case some columns read as character, coerce to numeric safely
for (cc in final_coverm_cols) {
  feat[[cc]] <- suppressWarnings(as.numeric(feat[[cc]]))
}
feat[is.na(feat)] <- 0

# Prepare output matrix:
# row1: class + classes
# row2: SampleID + sample IDs
# row3+: Genome + abundances
out_mat <- rbind(
  c("class", final_class),
  c("SampleID", final_sample_ids),
  cbind(as.character(feat$Genome), as.matrix(feat[, final_coverm_cols, drop = FALSE]))
)

out_tsv <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/lefse_input_repMAGs_tpm.tsv"

# write (tab-delimited, no quotes)
write.table(
  out_mat,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("✅ LEfSe 输入文件已生成：", out_tsv)
message("   使用的丰度列：", abund_suffix)
message("   样本数：", length(final_sample_ids), "；MAG 数：", nrow(feat))

otu_tpm <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/lefse_input_repMAGs_tpm.tsv", skip = 1)

diff_mag <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 8-AEE PEE differential MAGs.xlsx")
diff_mag <- diff_mag %>% subset(Species != 'Others')
diff_mag <- diff_mag$MAG_formal

diff_tpm <- coverm %>% subset(coverm$MAG_formal %in% diff_mag)

diff_tpm <- diff_tpm %>% select(contains('MAG_formal'), contains('TPM'))

openxlsx::write.xlsx(diff_tpm, "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/GWAS_MAG_TPM.xlsx")
