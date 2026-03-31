#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(colorspace)
  library(purrr)
})

# =========================
# 0) 文件路径
# =========================
metab_file <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/代谢组数据/归一化代谢组数据.xlsx"
group_file <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx"

outdir <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/代谢组数据/Figures/AEE_PEE_pathway_metabolites_violin"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) 读取分组信息
# =========================
group_df <- read_xlsx(group_file, sheet = "VFA_rumenotype_assignments") %>%
  transmute(
    SampleID   = as.character(SampleID),
    Rumenotype = case_when(
      Rumenotype %in% c("AEE", "PEE") ~ Rumenotype,
      Rumenotype %in% c("Acetate_type", "Acetate Type", "acetate_type") ~ "AEE",
      Rumenotype %in% c("Propionate_type", "Propionate Type", "propionate_type") ~ "PEE",
      TRUE ~ as.character(Rumenotype)
    )
  ) %>%
  filter(Rumenotype %in% c("AEE", "PEE")) %>%
  distinct()

# =========================
# 2) 读取代谢组数据
# =========================
metab_raw <- read_xlsx(metab_file, sheet = 1)

if (!"Metabolite" %in% colnames(metab_raw)) {
  stop("代谢组表中未找到 'Metabolite' 列，请检查表头。")
}

sample_cols <- intersect(colnames(metab_raw), group_df$SampleID)
if (length(sample_cols) == 0) {
  stop("代谢组表中没有找到与分组表 SampleID 对应的样本列。")
}

message("Matched sample columns: ", length(sample_cols))

# =========================
# 3) 定义通路代谢物匹配表
#    Pattern 用正则表达式匹配真实代谢物名
# =========================
target_tbl <- tibble::tribble(
  ~Branch,              ~Node,                         ~Target,                  ~Pattern,
  "Shared upstream",    "Fermentable carbohydrates",   "Sucrose",                "^Sucrose$",
  "Shared upstream",    "Fermentable carbohydrates",   "Cellobiose",             "^Cellobiose$",
  "Shared upstream",    "Fermentable carbohydrates",   "Maltotetraose",          "^Maltotetraose$",
  "Shared upstream",    "Fermentable carbohydrates",   "Melezitose",             "^Melezitose$",
  "Shared upstream",    "Fermentable carbohydrates",   "Verbascose",             "^Verbascose$",
  "Shared upstream",    "Hexose input",                "Glucose",                "^D-Glucose$|^Alpha-d-glucose$|^Glucose$",
  
  "Shared upstream",    "Glycolysis",                  "3-Phosphoglycerate",     "^3-Phosphoglycerate$",
  "Shared upstream",    "Glycolysis",                  "PEP",                    "^Phosphoenol pyruvate$|^Phosphoenolpyruvate$",
  "Shared upstream",    "Glycolysis",                  "Pyruvate",               "^Pyruvate$",
  
  "AEE branch",         "Acetate route",               "Acetyl-CoA",             "^Acetyl-CoA$",
  "AEE branch",         "Acetate route",               "Acetyl-phosphate",       "^Acetyl phosphate$|^Acetyl-phosphate$",
  "AEE branch",         "Acetate route",               "Acetate",                "^Acetate$|^Acetic acid$",
  
  "PEE branch",         "Succinate-propionate route",  "Oxaloacetate",           "^Oxaloacetate$",
  "PEE branch",         "Succinate-propionate route",  "Malate",                 "^L-malic acid$|^Malic acid$|^Malate$",
  "PEE branch",         "Succinate-propionate route",  "Fumarate",               "^Fumaric acid$|^Fumarate$",
  "PEE branch",         "Succinate-propionate route",  "Succinate",              "^Succinic acid$|^Succinate$",
  "PEE branch",         "Succinate-propionate route",  "Methylmalonate",         "^Methylmalonic Acid$",
  "PEE branch",         "Succinate-propionate route",  "Propionate",             "^Propionic Acid$",
  
  "Amino acid input",   "Asp family",                  "Aspartate",              "^L-Aspartic Acid$|^N-Acetyl-L-aspartic acid$",
  "Amino acid input",   "Homoserine / Thr family",     "Homoserine",             "^DL-Homoserine$|^L-Homoserine$|^O-Acetyl-L-homoserine$",
  "Amino acid input",   "Homoserine / Thr family",     "N-acetylthreonine",      "^N-acetylthreonine$",
  "Amino acid input",   "Homoserine / Thr family",     "2-aminobutyric acid",    "^2-aminobutyric acid$",
  "Amino acid input",   "Ile / Val family",            "Isoleucine",             "^L-Isoleucine$|^N-Acetyl-D-alloisoleucine$",
  "Amino acid input",   "Ile / Val family",            "Valine",                 "^L-Valine$",
  "Amino acid input",   "Other amino acids",           "Glycine",                "^L-Glycine$",
  "Amino acid input",   "Other amino acids",           "Proline",                "^L-Proline$",
  
  "Amino acid input",   "Asp/Glu-containing peptides", "Asp-Glu",                "^Asp-Glu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Asp-Thr",                "^Asp-Thr$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Asp-Glu-Leu",            "^Asp-Glu-Leu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Asp-Asn",                "^Asp-Asn$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Val-Glu-Asp",            "^Val-Glu-Asp$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Asn-Glu-Asp",            "^Asn-Glu-Asp$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Asp-Leu",            "^Glu-Asp-Leu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Leu-Glu-Asp",            "^Leu-Glu-Asp$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Ile-Asp-Asp",            "^Ile-Asp-Asp$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Gly-Glu",                "^Gly-Glu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Thr-Glu",                "^Thr-Glu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Glu",                "^Glu-Glu$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Pro-Lys",            "^Glu-Pro-Lys$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Phe-Asp",            "^Glu-Phe-Asp$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Glu-Asn",            "^Glu-Glu-Asn$",
  "Amino acid input",   "Asp/Glu-containing peptides", "Glu-Glu-Phe",            "^Glu-Glu-Phe$"
)

target_levels <- unique(target_tbl$Target)

# =========================
# 4) 在代谢组表中匹配真实代谢物名
# =========================
all_metabolites <- unique(as.character(metab_raw$Metabolite))

match_tbl <- target_tbl %>%
  rowwise() %>%
  mutate(
    MatchedMetabolites = list(all_metabolites[str_detect(all_metabolites, regex(Pattern, ignore_case = TRUE))]),
    n_match = length(MatchedMetabolites),
    MatchedText = ifelse(n_match == 0, NA_character_, paste(MatchedMetabolites, collapse = " | "))
  ) %>%
  ungroup()

matched_tbl <- match_tbl %>% filter(n_match > 0)
missing_tbl <- match_tbl %>% filter(n_match == 0)

write.csv(matched_tbl %>% select(Branch, Node, Target, Pattern, n_match, MatchedText),
          file = file.path(outdir, "matched_targets.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(missing_tbl %>% select(Branch, Node, Target, Pattern),
          file = file.path(outdir, "missing_targets.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("Matched targets: ", nrow(matched_tbl))
message("Missing targets: ", nrow(missing_tbl))

if (nrow(matched_tbl) == 0) {
  stop("一个目标代谢物都没有匹配到，请检查代谢物名称。")
}

plot_map <- matched_tbl %>%
  select(Branch, Node, Target, MatchedMetabolites) %>%
  tidyr::unnest(MatchedMetabolites) %>%
  rename(Metabolite = MatchedMetabolites)

# =========================
# 5) 宽表转长表，并合并分组
# =========================
plot_df <- metab_raw %>%
  select(Metabolite, all_of(sample_cols)) %>%
  # filter(Metabolite %in% plot_map$Metabolite) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  left_join(group_df, by = "SampleID") %>%
  filter(Rumenotype %in% c("AEE", "PEE")) %>%
  left_join(plot_map, by = "Metabolite") %>%
  group_by(Branch, Node, Target, Metabolite, SampleID, Rumenotype) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Rumenotype = factor(Rumenotype, levels = c("AEE", "PEE")),
    Target = factor(Target, levels = target_levels),
    facet_lab = ifelse(Target == Metabolite,
                       as.character(Target),
                       paste0(as.character(Target), "\n[", Metabolite, "]"))
  )

# =========================
# 6) 统计检验（每个代谢物 Wilcoxon）
# =========================
safe_wilcox <- function(df) {
  out <- tryCatch({
    wilcox.test(Abundance ~ Rumenotype, data = df)$p.value
  }, error = function(e) NA_real_)
  out
}

stat_df <- plot_df %>%
  group_by(Branch, Node, Target, Metabolite, facet_lab) %>%
  summarise(
    p_value = safe_wilcox(cur_data()),
    ymax = max(Abundance, na.rm = TRUE),
    ymin = min(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    yrange = ifelse(is.finite(ymax - ymin) & (ymax - ymin) > 0, ymax - ymin, abs(ymax) * 0.1 + 0.1),
    y_pos = ymax + 0.12 * yrange,
    p_label = case_when(
      is.na(p_value)      ~ "NA",
      p_value < 0.001     ~ "***",
      p_value < 0.01      ~ "**",
      p_value < 0.05      ~ "*",
      TRUE                ~ paste0("P = ", signif(p_value, 2))
    )
  )

write.csv(stat_df,
          file = file.path(outdir, "AEE_PEE_metabolite_wilcox_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# =========================
# 7) 作图函数
# =========================
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA)
  )

safe_name <- function(x) {
  x %>%
    str_replace_all("[/\\\\:*?\"<>|]", "_") %>%
    str_replace_all("\\s+", "_")
}
# =========================
# 10) 每个代谢物单独输出一张图
# =========================
single_dir <- file.path(outdir, "single_metabolite_plots")
dir.create(single_dir, showWarnings = FALSE)

single_meta_list <- split(plot_df, plot_df$Metabolite)

nm <- "L-Valine, N-(2-hydroxy-3-butenyl)-"

  d <- single_meta_list[[nm]]
  
  p <- ggplot(d, aes(x = Rumenotype, y = Abundance, fill = Rumenotype, color = Rumenotype)) +
    # geom_violin(trim = FALSE, alpha = 0.35, linewidth = 0.45, width = 0.4, na.rm = TRUE) +
    # geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.9, linewidth = 0.5, na.rm = TRUE) +
    # 均值水平线段
    stat_summary(
      fun = mean,
      fun.min = median,
      fun.max = median,
      geom = "bar",
      width = 0.3,       # 控制线段横向延伸长度
      linewidth = 0.5,        # 线段粗细
      show.legend = FALSE,
      position = position_dodge(width = 0.1)
    )+
    geom_jitter(width = 0.1, size = 1.2, alpha = 0.8, na.rm = TRUE) +
    scale_fill_manual(values = c("AEE" = "#2e528f",
                                 "PEE" = "#b41d23")) +
    scale_color_manual(values = c("AEE" = darken("#2e528f", 0.5),
                                  "PEE" = darken("#b41d23", 0.5))) +
    # base_theme +
    theme_void()+
    theme(
      legend.position = "none",
    )
    
    
  p
  
  fn <- safe_name(nm)
  # fn
  ggsave(file.path(single_dir, paste0(fn, ".pdf")), p,
         width = 3, height = 2, units = "in", dpi = 300)


message("Done.")
message("Output directory: ", outdir)
  
  
  
VFA <- read_xlsx(group_file, sheet = "VFA_proportion")

ggplot(VFA, aes(x = Rumenotype, y = Propionic_acid, fill = Rumenotype, color = Rumenotype)) +
  # geom_violin(trim = FALSE, alpha = 0.35, linewidth = 0.45, width = 0.4, na.rm = TRUE) +
  # geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.9, linewidth = 0.5, na.rm = TRUE) +
  # 均值水平线段
  stat_summary(
    fun = mean,
    fun.min = median,
    fun.max = median,
    geom = "bar",
    width = 0.3,       # 控制线段横向延伸长度
    linewidth = 0.5,        # 线段粗细
    show.legend = FALSE,
    position = position_dodge(width = 0.1)
  )+
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.8, na.rm = TRUE) +
  scale_fill_manual(values = c("AEE" = "#2e528f",
                               "PEE" = "#b41d23")) +
  scale_color_manual(values = c("AEE" = darken("#2e528f", 0.5),
                                "PEE" = darken("#b41d23", 0.5))) +
  # base_theme +
  theme_void()+
  theme(
    legend.position = "none",
  )#+
  # ggsignif::geom_signif(comparisons = list(c('AEE','PEE')))

ggsave(file.path(single_dir, paste0("Propionic_acid", ".pdf")),
       width = 3, height = 2, units = "in", dpi = 300)



assign <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx")

EC <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/EC/TPM_by_EC.tsv")

KO <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/EC/TPM_by_KEGG_ko.tsv")

# =========================
# 在你现有代码后面续写
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(openxlsx)
})

# -------------------------
# 1. 规范分组表
# -------------------------
assign2 <- assign

# 如果你前面没有指定sheet，这里强制重新读一下指定sheet更稳
if (!("Rumenotype" %in% names(assign2))) {
  assign2 <- openxlsx::read.xlsx(
    "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx",
    sheet = "VFA_rumenotype_assignments"
  )
}

# 自动识别样本列和分组列
sample_col <- c("SampleID", "sampleID", "Sample", "sample", "ID", "id")
group_col  <- c("Rumenotype", "rumenotype", "Group", "group", "Type", "type")

sample_col <- sample_col[sample_col %in% names(assign2)][1]
group_col  <- group_col[group_col %in% names(assign2)][1]

if (is.na(sample_col) || is.na(group_col)) {
  stop("assign 表中未识别到 SampleID 或 Rumenotype 列，请检查列名。")
}

group_df <- assign2 %>%
  transmute(
    SampleID = as.character(.data[[sample_col]]),
    Rumenotype = as.character(.data[[group_col]])
  ) %>%
  mutate(
    Rumenotype = case_when(
      Rumenotype %in% c("AEE", "PEE") ~ Rumenotype,
      str_detect(Rumenotype, regex("acetate", ignore_case = TRUE)) ~ "AEE",
      str_detect(Rumenotype, regex("propionate", ignore_case = TRUE)) ~ "PEE",
      TRUE ~ Rumenotype
    )
  ) %>%
  filter(Rumenotype %in% c("AEE", "PEE")) %>%
  distinct()

cat("Group table loaded: ", nrow(group_df), " samples\n")

# -------------------------
# 2. 我前面列出的目标 EC 和 KO
#    优先保留主酶和常见替代酶
# -------------------------

target_EC <- tibble::tribble(
  ~Pathway_step, ~Enzyme, ~EC_ID,
  "Sucrose -> glucose/fructose", "beta-fructofuranosidase / invertase", "3.2.1.26",
  "Cellobiose -> glucose", "beta-glucosidase", "3.2.1.21",
  "Maltotetraose -> smaller sugars", "alpha-amylase", "3.2.1.1",
  "Maltotetraose -> glucose", "alpha-glucosidase", "3.2.1.20",
  "Melezitose/Verbascose degradation", "alpha-galactosidase", "3.2.1.22",
  
  "Glucose -> G6P", "hexokinase", "2.7.1.1",
  "Glucose -> G6P", "glucokinase", "2.7.1.2",
  "G6P <-> F6P", "glucose-6-phosphate isomerase", "5.3.1.9",
  "F6P -> F1,6BP", "6-phosphofructokinase", "2.7.1.11",
  "F1,6BP -> GAP + DHAP", "fructose-bisphosphate aldolase", "4.1.2.13",
  "DHAP <-> GAP", "triosephosphate isomerase", "5.3.1.1",
  "GAP -> 1,3-BPG", "glyceraldehyde-3-phosphate dehydrogenase", "1.2.1.12",
  "GAP -> 1,3-BPG", "glyceraldehyde-3-phosphate dehydrogenase", "1.2.1.59",
  "1,3-BPG -> 3-PG", "phosphoglycerate kinase", "2.7.2.3",
  "3-PG <-> 2-PG", "phosphoglycerate mutase", "5.4.2.11",
  "3-PG <-> 2-PG", "phosphoglycerate mutase", "5.4.2.12",
  "2-PG -> PEP", "enolase", "4.2.1.11",
  "PEP -> pyruvate", "pyruvate kinase", "2.7.1.40",
  
  "Pyruvate -> acetyl-CoA", "pyruvate:ferredoxin oxidoreductase", "1.2.7.1",
  "Acetyl-CoA -> acetyl-phosphate", "phosphate acetyltransferase", "2.3.1.8",
  "Acetyl-phosphate -> acetate", "acetate kinase", "2.7.2.1",
  
  "PEP -> oxaloacetate", "phosphoenolpyruvate carboxylase", "4.1.1.31",
  "PEP -> oxaloacetate", "phosphoenolpyruvate carboxykinase", "4.1.1.49",
  "PEP -> oxaloacetate", "phosphoenolpyruvate carboxykinase", "4.1.1.32",
  "Pyruvate -> oxaloacetate", "pyruvate carboxylase", "6.4.1.1",
  "Aspartate <-> oxaloacetate", "aspartate aminotransferase", "2.6.1.1",
  "Oxaloacetate -> malate", "malate dehydrogenase", "1.1.1.37",
  "Oxaloacetate -> malate", "malate dehydrogenase (quinone)", "1.1.5.4",
  "Malate <-> fumarate", "fumarase", "4.2.1.2",
  "Fumarate -> succinate", "fumarate reductase / succinate dehydrogenase", "1.3.5.4",
  "Succinate <-> succinyl-CoA", "succinyl-CoA synthetase", "6.2.1.5",
  "Succinyl-CoA <-> methylmalonyl-CoA", "methylmalonyl-CoA mutase", "5.4.99.2",
  "Methylmalonyl-CoA epimerization", "methylmalonyl-CoA epimerase", "5.1.99.1",
  "Methylmalonyl-CoA decarboxylation", "methylmalonyl-CoA decarboxylase", "4.1.1.41",
  "Propionyl-CoA -> propionate", "propionate CoA-transferase", "2.8.3.1",
  "Propionyl-phosphate -> propionate", "propionate kinase", "2.7.2.15",
  
  "Peptide hydrolysis", "peptidases", "3.4.-.-",
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "1.4.1.2",
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "1.4.1.3",
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "1.4.1.4",
  "Threonine -> 2-ketobutyrate", "threonine dehydratase", "4.3.1.19",
  "Val/Ile transamination", "branched-chain amino acid aminotransferase", "2.6.1.42",
  "Branched-chain 2-oxoacid -> acyl-CoA", "branched-chain 2-oxoacid dehydrogenase", "1.2.4.4",
  "Branched-chain 2-oxoacid -> acyl-CoA", "dihydrolipoamide acyltransferase", "2.3.1.168",
  
  "OAA + acetyl-CoA -> citrate", "citrate synthase", "2.3.3.1",
  "Citrate <-> isocitrate", "aconitate hydratase", "4.2.1.3",
  "Isocitrate -> 2-oxoglutarate", "isocitrate dehydrogenase", "1.1.1.42",
  "Isocitrate -> 2-oxoglutarate", "isocitrate dehydrogenase", "1.1.1.41",
  "2-oxoglutarate -> succinyl-CoA", "2-oxoglutarate dehydrogenase", "1.2.4.2",
  "2-oxoglutarate -> succinyl-CoA", "dihydrolipoamide succinyltransferase", "2.3.1.61",
  "2-oxoglutarate -> succinyl-CoA", "dihydrolipoamide dehydrogenase", "1.8.1.4"
)

target_KO <- tibble::tribble(
  ~Pathway_step, ~Enzyme, ~KO_ID,
  "Sucrose -> glucose/fructose", "beta-fructofuranosidase / invertase", "K01193",
  "Cellobiose -> glucose", "beta-glucosidase", "K01188",
  "Cellobiose -> glucose", "beta-glucosidase", "K05349",
  "Cellobiose -> glucose", "beta-glucosidase", "K05350",
  "Maltotetraose -> smaller sugars", "alpha-amylase", "K01176",
  "Maltotetraose -> glucose", "alpha-glucosidase", "K01187",
  "Melezitose/Verbascose degradation", "alpha-galactosidase", "K07406",
  
  "Glucose -> G6P", "hexokinase", "K00844",
  "Glucose -> G6P", "glucokinase", "K12407",
  "Glucose -> G6P", "glucokinase", "K00845",
  "G6P <-> F6P", "glucose-6-phosphate isomerase", "K01810",
  "F6P -> F1,6BP", "6-phosphofructokinase", "K00850",
  "F1,6BP -> GAP + DHAP", "fructose-bisphosphate aldolase", "K01623",
  "F1,6BP -> GAP + DHAP", "fructose-bisphosphate aldolase", "K01624",
  "DHAP <-> GAP", "triosephosphate isomerase", "K01803",
  "GAP -> 1,3-BPG", "glyceraldehyde-3-phosphate dehydrogenase", "K00134",
  "GAP -> 1,3-BPG", "glyceraldehyde-3-phosphate dehydrogenase", "K00150",
  "1,3-BPG -> 3-PG", "phosphoglycerate kinase", "K00927",
  "3-PG <-> 2-PG", "phosphoglycerate mutase", "K01834",
  "3-PG <-> 2-PG", "phosphoglycerate mutase", "K15634",
  "3-PG <-> 2-PG", "phosphoglycerate mutase", "K15633",
  "2-PG -> PEP", "enolase", "K01689",
  "PEP -> pyruvate", "pyruvate kinase", "K00873",
  
  "Pyruvate -> acetyl-CoA", "pyruvate:ferredoxin oxidoreductase", "K00169",
  "Pyruvate -> acetyl-CoA", "pyruvate:ferredoxin oxidoreductase", "K00170",
  "Pyruvate -> acetyl-CoA", "pyruvate:ferredoxin oxidoreductase", "K00171",
  "Pyruvate -> acetyl-CoA", "pyruvate:ferredoxin oxidoreductase", "K03737",
  "Acetyl-CoA -> acetyl-phosphate", "phosphate acetyltransferase", "K00625",
  "Acetyl-CoA -> acetyl-phosphate", "phosphate acetyltransferase", "K13788",
  "Acetyl-CoA -> acetyl-phosphate", "phosphate acetyltransferase", "K15024",
  "Acetyl-phosphate -> acetate", "acetate kinase", "K00925",
  
  "PEP -> oxaloacetate", "phosphoenolpyruvate carboxylase", "K01595",
  "PEP -> oxaloacetate", "ATP-dependent phosphoenolpyruvate carboxykinase", "K01610",
  "PEP -> oxaloacetate", "GTP-dependent phosphoenolpyruvate carboxykinase", "K01596",
  "Pyruvate -> oxaloacetate", "pyruvate carboxylase", "K01958",
  "Pyruvate -> oxaloacetate", "pyruvate carboxylase", "K01959",
  "Pyruvate -> oxaloacetate", "pyruvate carboxylase", "K01960",
  "Aspartate <-> oxaloacetate", "aspartate aminotransferase", "K00812",
  "Aspartate <-> oxaloacetate", "aspartate aminotransferase", "K00813",
  "Oxaloacetate -> malate", "malate dehydrogenase", "K00024",
  "Oxaloacetate -> malate", "malate dehydrogenase", "K00025",
  "Oxaloacetate -> malate", "malate dehydrogenase (quinone)", "K00116",
  "Malate <-> fumarate", "fumarase", "K01676",
  "Fumarate -> succinate", "fumarate reductase subunit A", "K00244",
  "Fumarate -> succinate", "fumarate reductase subunit B", "K00245",
  "Fumarate -> succinate", "fumarate reductase subunit C", "K00246",
  "Fumarate -> succinate", "fumarate reductase subunit D", "K00247",
  "Succinate <-> succinyl-CoA", "succinyl-CoA synthetase beta", "K01902",
  "Succinate <-> succinyl-CoA", "succinyl-CoA synthetase alpha", "K01903",
  "Succinyl-CoA <-> methylmalonyl-CoA", "methylmalonyl-CoA mutase", "K01847",
  "Succinyl-CoA <-> methylmalonyl-CoA", "methylmalonyl-CoA mutase", "K01848",
  "Succinyl-CoA <-> methylmalonyl-CoA", "methylmalonyl-CoA mutase", "K01849",
  "Methylmalonyl-CoA decarboxylation", "methylmalonyl-CoA decarboxylase", "K11264",
  "Methylmalonyl-CoA carboxyltransferase", "methylmalonyl-CoA carboxyltransferase", "K17489",
  "Methylmalonyl-CoA carboxyltransferase", "methylmalonyl-CoA carboxyltransferase", "K17490",
  "Propionyl-CoA -> propionate", "propionate CoA-transferase", "K01026",
  "Propionyl-phosphate -> propionate", "phosphate propanoyltransferase", "K13923",
  "Propionyl-phosphate -> propionate", "propionate kinase", "K19697",
  "Propionyl-phosphate -> propionate", "propionate kinase", "K00932",
  
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "K00260",
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "K00261",
  "Glutamate -> 2-oxoglutarate", "glutamate dehydrogenase", "K00262",
  "Threonine -> 2-ketobutyrate", "threonine dehydratase", "K01754",
  "Val/Ile transamination", "branched-chain amino acid aminotransferase", "K00826",
  "Branched-chain 2-oxoacid -> acyl-CoA", "branched-chain 2-oxoacid dehydrogenase", "K00166",
  "Branched-chain 2-oxoacid -> acyl-CoA", "branched-chain 2-oxoacid dehydrogenase", "K00167",
  "Branched-chain 2-oxoacid -> acyl-CoA", "branched-chain 2-oxoacid dehydrogenase", "K09699",
  
  "OAA + acetyl-CoA -> citrate", "citrate synthase", "K01647",
  "Citrate <-> isocitrate", "aconitate hydratase", "K01681",
  "Citrate <-> isocitrate", "aconitate hydratase", "K01682",
  "Isocitrate -> 2-oxoglutarate", "isocitrate dehydrogenase", "K00031",
  "Isocitrate -> 2-oxoglutarate", "isocitrate dehydrogenase", "K00030",
  "2-oxoglutarate -> succinyl-CoA", "2-oxoglutarate dehydrogenase E1", "K00164",
  "2-oxoglutarate -> succinyl-CoA", "dihydrolipoamide succinyltransferase E2", "K00658",
  "2-oxoglutarate -> succinyl-CoA", "dihydrolipoamide dehydrogenase E3", "K00382"
)

# -------------------------
# 3. 编号规范化函数
# -------------------------
norm_ec <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_remove("^ec:") %>%
    str_remove("^EC:") %>%
    str_replace_all("\\s+", "")
}

norm_ko <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_remove("^ko:") %>%
    str_remove("^KO:") %>%
    str_to_upper() %>%
    str_extract("K\\d{5}")
}

# -------------------------
# 4. TPM表转成长表
#    自动判断是 feature × sample 还是 sample × feature
# -------------------------
make_long_tpm <- function(df, group_df, type = c("EC", "KO")) {
  type <- match.arg(type)
  df <- as.data.frame(df)
  
  # 样本ID
  sample_ids <- unique(group_df$SampleID)
  
  # 判断列名里是否有样本
  sample_cols_in_colnames <- intersect(colnames(df), sample_ids)
  
  if (length(sample_cols_in_colnames) >= 2) {
    # 说明是 feature × sample
    id_candidates <- switch(
      type,
      EC = c("EC", "ec", "Feature", "feature", "ID", "id", colnames(df)[1]),
      KO = c("KO", "ko", "KEGG_ko", "KEGG_KO", "Feature", "feature", "ID", "id", colnames(df)[1])
    )
    id_col <- id_candidates[id_candidates %in% colnames(df)][1]
    
    out <- df %>%
      dplyr::select(all_of(c(id_col, sample_cols_in_colnames))) %>%
      dplyr::rename(feature_raw = all_of(id_col)) %>%
      tidyr::pivot_longer(
        cols = -feature_raw,
        names_to = "SampleID",
        values_to = "TPM"
      )
  } else {
    # 说明可能是 sample × feature
    sample_candidates <- c("SampleID", "sampleID", "Sample", "sample", "ID", "id", colnames(df)[1])
    sample_col2 <- sample_candidates[sample_candidates %in% colnames(df)][1]
    
    out <- df %>%
      dplyr::rename(SampleID = all_of(sample_col2)) %>%
      dplyr::filter(SampleID %in% sample_ids) %>%
      tidyr::pivot_longer(
        cols = -SampleID,
        names_to = "feature_raw",
        values_to = "TPM"
      )
  }
  
  out <- out %>%
    mutate(
      feature_raw = as.character(feature_raw),
      TPM = as.numeric(TPM),
      feature_id = if (type == "EC") norm_ec(feature_raw) else norm_ko(feature_raw)
    ) %>%
    left_join(group_df, by = "SampleID") %>%
    filter(!is.na(Rumenotype))
  
  out
}

EC_long <- make_long_tpm(EC, group_df, type = "EC")
KO_long <- make_long_tpm(KO, group_df, type = "KO")

cat("EC_long rows:", nrow(EC_long), "\n")
cat("KO_long rows:", nrow(KO_long), "\n")

# -------------------------
# 5. Wilcoxon差异分析函数
# -------------------------
run_wilcox_target <- function(long_df, target_df, id_col = c("EC_ID", "KO_ID")) {
  id_col <- match.arg(id_col)
  
  target_df2 <- target_df %>%
    mutate(
      feature_id = if (id_col == "EC_ID") norm_ec(.data[[id_col]]) else norm_ko(.data[[id_col]])
    ) %>%
    distinct()
  
  # 只保留目标编号
  dat <- long_df %>%
    inner_join(target_df2 %>% select(feature_id), by = "feature_id")
  
  # 每个 feature 做 wilcox
  stat_df <- dat %>%
    group_by(feature_id) %>%
    summarise(
      n_AEE = sum(Rumenotype == "AEE" & !is.na(TPM)),
      n_PEE = sum(Rumenotype == "PEE" & !is.na(TPM)),
      mean_AEE = mean(TPM[Rumenotype == "AEE"], na.rm = TRUE),
      mean_PEE = mean(TPM[Rumenotype == "PEE"], na.rm = TRUE),
      median_AEE = median(TPM[Rumenotype == "AEE"], na.rm = TRUE),
      median_PEE = median(TPM[Rumenotype == "PEE"], na.rm = TRUE),
      wilcox_p = tryCatch(
        wilcox.test(TPM ~ Rumenotype)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      padj_BH = p.adjust(wilcox_p, method = "BH"),
      direction = case_when(
        is.na(median_AEE) | is.na(median_PEE) ~ NA_character_,
        median_AEE > median_PEE ~ "AEE_higher",
        median_AEE < median_PEE ~ "PEE_higher",
        TRUE ~ "equal"
      ),
      log2FC_AEE_vs_PEE = log2((mean_AEE + 1e-8) / (mean_PEE + 1e-8))
    )
  
  # 合并回注释表
  res <- target_df2 %>%
    left_join(stat_df, by = "feature_id") %>%
    arrange(padj_BH, wilcox_p)
  
  res
}

EC_res <- run_wilcox_target(EC_long, target_EC, id_col = "EC_ID")
KO_res <- run_wilcox_target(KO_long, target_KO, id_col = "KO_ID")

# -------------------------
# 6. 输出结果
# -------------------------
outdir <- "D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/AEE_PEE_target_enzyme_wilcox"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(EC_res, file.path(outdir, "EC_wilcox_results.csv"))
readr::write_csv(KO_res, file.path(outdir, "KO_wilcox_results.csv"))

openxlsx::write.xlsx(
  list(
    Group_info = group_df,
    Target_EC_list = target_EC,
    Target_KO_list = target_KO,
    EC_wilcox_results = EC_res,
    KO_wilcox_results = KO_res
  ),
  file = file.path(outdir, "AEE_PEE_target_enzyme_wilcox_results.xlsx"),
  overwrite = TRUE
)

# -------------------------
# 7. 简单打印显著结果
# -------------------------
cat("\n================ EC significant (BH < 0.05) ================\n")
print(EC_res %>% filter(!is.na(padj_BH), padj_BH < 0.05))

cat("\n================ KO significant (BH < 0.05) ================\n")
print(KO_res %>% filter(!is.na(padj_BH), padj_BH < 0.05))

cat("\nResults saved to:\n", outdir, "\n")



#########################
# 相关性结果
#########################

KO_AP <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 9-KO association with AP.xlsx")

KO_rumenotype <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 10-KO association with Rumenotype.xlsx")

KO_AP_sig <- KO_AP %>% subset(pvalue_AP < 0.05)

KO_rumenotype_sig <- KO_rumenotype %>% subset(pvalue_group < 0.05)

KO_sig <- full_join(KO_AP, KO_rumenotype, by = 'KO')

KO_sig$consistent <- KO_sig$beta_AP_CLR*KO_sig$beta_group_CLR 

KO_sig$consistentency <- ifelse(KO_sig$consistent >0, 'consistent', 'inconsistent')

KO_sig <- KO_sig %>% subset(consistentency == 'consistent')

KO_sig <- KO_sig %>% subset(pvalue_AP < 0.05 | pvalue_group < 0.05)

KO_sig$KO <- sub("^ko:", "", KO_sig$KO)

table(KO_sig$KO %in% target_KO$KO_ID)

