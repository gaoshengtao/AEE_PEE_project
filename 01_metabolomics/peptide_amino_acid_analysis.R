## =========================================================
##  Peptide/AA extraction + length stats + density plot + AA composition
## =========================================================

## 0) Packages -------------------------------------------------------------
pkgs <- c("readxl","dplyr","tidyr","stringr","purrr","ggplot2","forcats")
invisible(lapply(pkgs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))
library(readxl); library(dplyr); library(tidyr); library(stringr)
library(purrr); library(ggplot2); library(forcats)

if (!requireNamespace("ggridges", quietly = TRUE)) install.packages("ggridges")
library(ggridges)

## 1) Read data ------------------------------------------------------------
infile <- "./FinalFigure/Supplemental file 2-AEE PEE differential metabolites.xlsx"

outdir <- "AEE_PEE_metabolomics_OPLSDA_Wilcox"

## 读第1个sheet；如果你的表在别的sheet，把 sheet= 改一下
df0 <- readxl::read_xlsx(infile, sheet = 1)

df0 <- df0 %>% subset(fermentation_class == 'Peptide- and amino acid–related fermentation products')

## 2) 自动识别列：代谢物名列、样本列、分组信息 --------------------------------
## 2.1 代谢物名列（尽量自动找：常见列名）
name_candidates <- c("Metabolite","metabolite","Compound","compound","Name","name",
                     "Metabolites","Metabolite_Name","Feature","feature","ID","id")
name_col <- intersect(name_candidates, names(df0))[1]
if (is.na(name_col)) {
  ## fallback：默认第一列为代谢物名
  name_col <- names(df0)[1]
}

## 2.2 分组信息：优先从“样本信息表”里取；如果你的 df0 本身就有 Group 列可直接用
group_candidates <- c("Group","Dir","group","Type","type","Rumenotype","rumenotype")
group_col_in_df <- intersect(group_candidates, names(df0))[1]

## 2.3 样本列：数值列（大多数代谢组表都是：前面注释列，后面样本丰度列）
num_cols <- names(df0)[sapply(df0, is.numeric)]

if (length(num_cols) == 0) {
  stop("没识别到数值型样本列：请检查 Excel 是否为字符型数字，或样本列不在当前sheet。")
}

## 如果你的表里还有一些数值注释列（VIP/logFC/pvalue等），可在这里手动排除
## 例如：
## num_cols <- setdiff(num_cols, c("VIP","logFC","pvalue"))

## 3) 氨基酸/肽段解析函数 ---------------------------------------------------
## 3.1 标准氨基酸三字母缩写集合 + 扩展一些常见
AA3 <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
         "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val",
         "Hyp","Cit","Orn","Abu","Nle","Nva")  # Hyp: hydroxyproline, Cit: citrulline, Orn: ornithine, Abu, norleucine(Nle), norvaline(Nva)

## 3.2 全名/变体 -> 3字母映射（可按你数据再补）
full_to_aa3 <- c(
  
  ## ====== 20种标准氨基酸（全名） ======
  "Glycine" = "Gly",
  "Alanine" = "Ala",
  "Valine" = "Val",
  "Leucine" = "Leu",
  "Isoleucine" = "Ile",
  "Serine" = "Ser",
  "Threonine" = "Thr",
  "Cysteine" = "Cys",
  "Methionine" = "Met",
  "Aspartic Acid" = "Asp",
  "Glutamic Acid" = "Glu",
  "Glutamic acid" = "Glu",
  "Asparagine" = "Asn",
  "Glutamine" = "Gln",
  "Lysine" = "Lys",
  "Arginine" = "Arg",
  "Histidine" = "His",
  "Phenylalanine" = "Phe",
  "Tyrosine" = "Tyr",
  "Tryptophan" = "Trp",
  "Proline" = "Pro",
  
  ## ====== D / L / DL 形式 ======
  "L-Glycine" = "Gly",
  "L-Alanine" = "Ala",
  "L-Valine" = "Val",
  "L-Leucine" = "Leu",
  "L-Isoleucine" = "Ile",
  "L-Serine" = "Ser",
  "L-Threonine" = "Thr",
  "L-Cysteine" = "Cys",
  "L-Methionine" = "Met",
  "L-Aspartic Acid" = "Asp",
  "L-Glutamic Acid" = "Glu",
  "L-Asparagine" = "Asn",
  "L-Glutamine" = "Gln",
  "L-Lysine" = "Lys",
  "L-Arginine" = "Arg",
  "L-Histidine" = "His",
  "L-Phenylalanine" = "Phe",
  "L-Tyrosine" = "Tyr",
  "L-Tryptophan" = "Trp",
  "L-Proline" = "Pro",
  
  "DL-Homoserine" = "Hse",
  "Homoserine" = "Hse",
  
  ## ====== 非标准氨基酸（你数据里出现的） ======
  "Norvaline" = "Nva",
  "L-norleucine" = "Nle",
  "Norleucine" = "Nle",
  "Ornithine" = "Orn",
  "Citrulline" = "Cit",
  "Hydroxyproline" = "Hyp",
  "Cis-4-Hydroxy-D-proline" = "Hyp",
  "Allysine" = "Aly",
  "Baikiain" = "Bai",
  "2-aminobutyric acid" = "Abu",
  "2-Aminobutyric Acid" = "Abu",
  
  ## ====== 常见修饰型 ======
  "N-acetylornithine" = "Orn",
  "N-Acetyl-L-Alanine" = "Ala",
  "N-Acetylasparagine" = "Asn",
  "N-Acetylaspartic acid" = "Asp",
  "N-acetyl-dl-glutamic acid" = "Glu",
  "N-Acetyl-DL-Valine" = "Val",
  "N-Acetyltryptophan" = "Trp",
  "N-Acetylglycine" = "Gly",
  "Acetylglycine" = "Gly",
  "Acetylleucine" = "Leu",
  
  ## ====== 特殊写法 ======
  "G-Glu-Val" = "Gly",  # G = Gly
  "Afalanine" = "Ala",
  
  ## ====== 其他在你数据中出现的氨基酸类 ======
  "Proline betaine" = "Pro",
  "Formiminoglutamic acid" = "Glu",
  "Ketoleucine" = "Leu",
  "4-Oxoproline" = "Pro",
  "6-Oxopiperidine-2-carboxylic acid" = "Lys",
  
  ## ====== 防止大小写问题 ======
  "glutamic acid" = "Glu",
  "aspartic acid" = "Asp"
)

## 3.3 一些“非连字符”的二肽/小肽名字（按你列表里常见的写了几个）
##     你如果有更多，可以继续加
special_pep <- c(
  
  ## ====== Digly / 拼写类 ======
  "Diglykokoll" = "Gly-Gly",
  
  ## ====== Histidyl- 系列 ======
  "Histidylglycine"     = "His-Gly",
  "Histidylthreonine"   = "His-Thr",
  "Histidyltyrosine"    = "His-Tyr",
  "Histidylvaline"      = "His-Val",
  
  ## ====== 其他 -yl 结构 ======
  "Aspartylglycosamine" = "Asp-Gly",
  "Arginine-glycine-aspartate-O-methyltyrosine amide" = "Arg-Gly-Asp-Tyr",
  
  ## ====== 你数据中可能被误识别的非标准书写 ======
  # "Tetrapeptide" = NA,   # 无法确定组成，不纳入
  
  ## ====== 修饰型但本质是二肽 ======
  "Sialorphin" = "Gln-His-Pro-Arg",   # 已知天然四肽
  
  ## ====== 其他明确是氨基酸缩合物 ======
  "Glycine, N-[1-(phenylacetyl)-L-prolyl]-" = "Pro-Gly",
  
  ## ====== 单字母开头缩写形式 ======
  "G-Glu-Val" = "Gly-Glu-Val"
)

## 3.4 清洗 + 提取肽段残基
extract_residues <- function(x){
  # 任何 NA / NULL / 空字符 / 非字符：直接返回空
  if (is.null(x) || length(x) == 0) return(character(0))
  x0 <- as.character(x)[1]
  if (is.na(x0)) return(character(0))
  x0 <- str_trim(x0)
  if (!nzchar(x0)) return(character(0))
  
  # special mapping
  if (exists("special_pep", inherits = TRUE) && x0 %in% names(special_pep)) {
    x0 <- special_pep[[x0]]
    if (is.na(x0) || !nzchar(x0)) return(character(0))
  }
  
  # 清洗
  x0 <- str_replace_all(x0, "\\s+", " ")
  x0 <- str_replace_all(x0, "^L-|^D-|^DL-|^N-|^O-", "")
  x0 <- str_replace_all(x0, "gamma-", "")
  x0 <- str_replace_all(x0, "epsilon-gamma-", "")
  x0 <- str_replace_all(x0, "Cyclo\\(|\\)", "")
  x0 <- str_replace_all(x0, "diketopiperazine", "")
  x0 <- str_replace_all(x0, "\\bAc\\b|N-Acetyl-|O-acetyl-|N-formyl-", "")
  x0 <- str_replace_all(x0, "[-_]+", "-")
  x0 <- str_trim(x0)
  
  # 全名 -> 3字母（如果你定义了 full_to_aa3）
  if (exists("full_to_aa3", inherits = TRUE)) {
    for (nm in names(full_to_aa3)) {
      x0 <- str_replace_all(x0, fixed(nm), full_to_aa3[[nm]])
    }
  }
  
  x0 <- str_trim(x0)
  if (!nzchar(x0)) return(character(0))
  
  # 解析 tokens
  has_dash <- isTRUE(str_detect(x0, "-"))
  toks <- if (has_dash) str_split(x0, "-", simplify = FALSE)[[1]] else x0
  toks <- toks[toks != ""]
  
  # 单字母简写
  toks <- ifelse(toks == "G", "Gly", toks)
  
  # 只保留AA3（如果你定义了 AA3）
  if (exists("AA3", inherits = TRUE)) {
    toks <- toks[toks %in% AA3]
  }
  
  toks
}

classify_compound <- function(name){
  name0 <- str_trim(as.character(name)[1])
  if (is.na(name0) || !nzchar(name0)) {
    return(list(type="Other", len=NA_integer_, residues=character(0)))
  }
  
  # 关键：Tetrapeptide 作为“4肽”，但残基未知
  if (name0 == "Tetrapeptide") {
    return(list(type="Peptide", len=4L, residues=character(0)))
  }
  
  aa <- extract_residues(name0)
  if (length(aa) == 1) return(list(type="AA", len=1L, residues=aa))
  if (length(aa) >= 2) return(list(type="Peptide", len=as.integer(length(aa)), residues=aa))
  list(type="Other", len=NA_integer_, residues=character(0))
}

## 4) 生成注释表（每个代谢物：类型/长度/残基） -------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(forcats)

## 0) 指定代谢物名列
name_col <- "Metabolite"

## 1) 生成注释表 annot（你前面已经解决过 classify_compound / extract_residues）
annot <- df0 %>%
  transmute(compound = .data[[name_col]]) %>%
  mutate(parsed = map(compound, classify_compound),
         type = map_chr(parsed, "type"),
         pep_len = map_int(parsed, "len"),
         residues = map(parsed, "residues")) %>%
  select(-parsed)

## 2) 合并到差异表
dat <- df0 %>%
  mutate(compound = .data[[name_col]]) %>%
  left_join(annot, by = "compound")

## 3) 只取 AA + Peptide
dat_ap <- dat %>%
  filter(type %in% c("AA","Peptide")) %>%
  mutate(
    len_bin = case_when(
      pep_len >= 8 ~ "8+",
      pep_len >= 1 ~ as.character(pep_len),
      TRUE ~ NA_character_
    ),
    len_bin = factor(len_bin, levels = c(as.character(1:7), "8+"))
  )

## =========================================================
## A) 不同长度肽段数量统计（差异列表层面）
## =========================================================
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggridges)

pep_only <- dat_ap %>%
  mutate(
    len_bin = ifelse(pep_len >= 8, "8+", as.character(pep_len)),
    len_num = ifelse(pep_len >= 8, 8, pep_len),
    # 规范方向：如果 Dir 已经是 AEE/PEE 就直接用；否则用 log2fc 判断
    Dir2 = case_when(
      Dir %in% c("AEE","PEE") ~ Dir,
      !is.na(log2fc) & log2fc > 0 ~ "AEE",
      !is.na(log2fc) & log2fc < 0 ~ "PEE",
      TRUE ~ NA_character_
    ),
    Dir2 = factor(Dir2, levels = c("AEE","PEE"))
  ) %>%
  filter(!is.na(Dir2))

len_count <- pep_only %>%
  count(Dir2, len_num, name = "n") %>%
  group_by(Dir2) %>%
  mutate(prop = n / sum(n)) %>%   # 组内相对比例
  ungroup()

p_line <- ggplot(len_count, 
                 aes(x = len_num, y = n, 
                     color = Dir2, group = Dir2)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = 1:8,
    labels = c("1","2","3","4","5","6","7","8+")
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = "Peptide length (number of amino acids)",
    y = "Number of peptides/amino acids",
    color = "Group"
  )+
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        legend.position = 'inside',
        legend.position.inside = c(0.8,0.8),
        # axis.text = element_blank()
        )+
  scale_color_manual(values = c('AEE'='#2e528f', 'PEE'='#b41d23'))

print(p_line)

ggsave(file.path(outdir, "肽段组成数量.pdf"), p_line, width = 6, height = 3)


## =========================================================
## B) 每个长度的小肽：氨基酸组成（Top10 + Others，不分组，count-based）
## =========================================================
if (!requireNamespace("colorspace", quietly = TRUE)) install.packages("colorspace")
library(colorspace)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1) 先把所有肽段的 residues 拆开，得到 AA 的总出现次数（用于选 Top10）
aa_all <- dat_ap %>%
  filter(type == "Peptide", !is.na(pep_len)) %>%
  unnest_longer(residues, values_to = "AA")

top10_aa <- aa_all %>%
  count(AA, sort = TRUE, name = "n") %>%
  slice_head(n = 10) %>%
  pull(AA)

# 2) 计算：每个长度的 AA 组成（Top10 + Others）
aa_comp_len_top10 <- dat_ap %>%
  mutate(
    len_bin = ifelse(pep_len >= 8, "8+", as.character(pep_len)),
    len_bin = factor(len_bin, levels = c(as.character(1:7), "8+"))
  ) %>%
  unnest_longer(residues, values_to = "AA") %>%
  mutate(
    AA_show = ifelse(AA %in% top10_aa, AA, "Others"),
    # 让图例顺序：Top10（按总体频率排序）+ Others
    AA_show = factor(
      AA_show,
      levels = c(top10_aa, "Others")
    )
  ) %>%
  count(len_bin, AA_show, name = "n") %>%
  group_by(len_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 3) 配色：Top10 用高级调色板（HCL），Others 灰色
pal_top10 <- colorspace::qualitative_hcl(length(top10_aa), palette = "Dark 3")
names(pal_top10) <- top10_aa
pal_all <- c(pal_top10, Others = "grey75")

p_aa_comp_top10 <- ggplot(aa_comp_len_top10, aes(x = len_bin, y = prop, fill = AA_show)) +
  geom_col(width = 0.85, color = NA) +
  scale_fill_manual(values = pal_all, drop = FALSE) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    legend.title = element_blank()
  ) +
  labs(
    x = "Peptide length (number of amino acids)",
    y = "Relative abundance (Top10 AA + Others)"
  )

print(p_aa_comp_top10)

ggsave(
  file.path(outdir, "不同长度肽段_氨基酸组成_Top10_Others.pdf"),
  p_aa_comp_top10, width = 7, height = 3
)


## =========================================================
## C) X=20AA(按总频次排序), Y=出现频次, 两组(AEE/PEE)并排堆叠(长度), 连续色阶
## =========================================================
if (!requireNamespace("viridisLite", quietly = TRUE)) install.packages("viridisLite")
library(viridisLite)

AA20 <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
          "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")

# 1) 生成：每个(Dir2, AA, 长度) 的出现次数
aa_len_grp <- dat_ap %>%
  mutate(
    len_num = ifelse(pep_len >= 8, 8L, as.integer(pep_len)),
    Dir2 = case_when(
      Dir %in% c("AEE","PEE") ~ Dir,
      !is.na(log2fc) & log2fc > 0 ~ "AEE",
      !is.na(log2fc) & log2fc < 0 ~ "PEE",
      TRUE ~ NA_character_
    ),
    Dir2 = factor(Dir2, levels = c("AEE","PEE"))
  ) %>%
  filter(!is.na(Dir2)) %>%
  unnest_longer(residues, values_to = "AA") %>%
  filter(AA %in% AA20) %>%
  count(AA, Dir2, len_num, name = "n")

# 2) 补全缺失组合为0（确保“有的长度某组没有=0”）
aa_len_grp <- aa_len_grp %>%
  tidyr::complete(
    AA = AA20,
    Dir2 = factor(c("AEE","PEE"), levels = c("AEE","PEE")),
    len_num = 1:8,
    fill = list(n = 0)
  )

# 3) 按总出现频次(AEE+PEE)排序 AA
aa_order <- aa_len_grp %>%
  group_by(AA) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(AA)

aa_len_grp <- aa_len_grp %>%
  mutate(
    AA = factor(AA, levels = aa_order),
    aa_idx = as.integer(AA),
    # 手动给 AEE/PEE 两根柱一个小偏移，实现“并排堆叠柱”
    x = aa_idx + ifelse(Dir2 == "AEE", -0.22, 0.22)
  )

# 4) 画图：每个氨基酸两根并排堆叠柱，堆叠色=长度(连续色阶)
p_aa_len_dodge_stack <- ggplot(aa_len_grp, aes(x = x, y = n, fill = len_num)) +
  geom_col(width = 0.38, color = NA) +
  scale_x_continuous(
    breaks = seq_along(aa_order),
    labels = aa_order,
    expand = expansion(add = 0.6)
  ) +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(8,'RdBu'),
    breaks = 1:8,
    labels = c("1","2","3","4","5","6","7","8+"),
    name = "Peptide length"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.length = unit(2, "mm")
  ) +
  labs(
    x = "Amino acids (ranked by total frequency in AEE+PEE)",
    y = "Occurrence count"
  )

print(p_aa_len_dodge_stack)

ggsave(file.path(outdir, "AA20_频次_按总频次排序_AEEvsPEE_长度堆叠_连续色阶.pdf"),
       p_aa_len_dodge_stack, width = 9, height = 3.5)
