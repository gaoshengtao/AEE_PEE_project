## =========================================================
## VFA profile 分型（Aitchison/CLR + PAM）
## 输入：VFA profiles.xlsx（SampleID + 6种VFA浓度）
## 输出：Rumenotype分型表、特征表、聚类评估图、PCoA散点图
## =========================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(cluster)
  library(clusterSim)
  library(ggplot2)
})

## -----------------------------
## 0) 读入数据
## -----------------------------
infile <- "./代谢组数据/VFA profiles.xlsx"  # 改成你的路径

raw <- read.xlsx(infile, sheet = 1)

# 期望的列名（按你表里的实际列名）
id_col <- "SampleID"
vfa_cols <- c("Acetic_acid","Propionic_acid","Isobutyric_acid",
              "Butyric_acid","Isovaleric_acid","Valeric_acid")


raw <- raw %>%
  rename_with(~trimws(.x)) %>%      # 清理列名空格（很重要）
  mutate(across(all_of(vfa_cols), as.numeric)) %>%
  group_by(.data[[id_col]]) %>%
  summarise(
    across(
      all_of(vfa_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )




## -----------------------------
## 1) 构建 VFA profile 特征
##    - prop: 6种VFA组成比例
##    - clr : compositional CLR 变换
##    - derived: TVFA, A:P, (A+B):P, BCVFA%
## -----------------------------
pseudo <- 1e-6

vfa_mat <- as.matrix(raw[, vfa_cols])
rownames(vfa_mat) <- raw[[id_col]]

tvfa <- rowSums(vfa_mat)
prop <- sweep(vfa_mat, 1, tvfa, "/")             # 组成比例

# CLR：对比例做clr（加pseudo避免0）
prop2 <- prop + pseudo
gm <- exp(rowMeans(log(prop2)))
clr <- log(prop2 / gm)

derived <- raw %>%
  transmute(
    SampleID = .data[[id_col]],
    TVFA = tvfa,
    AP   = (Acetic_acid + pseudo) / (Propionic_acid + pseudo),
    ABP  = (Acetic_acid + Butyric_acid + pseudo) / (Propionic_acid + pseudo),
    BCVFA = (Isobutyric_acid + Isovaleric_acid) / (TVFA + pseudo)
  )

## -----------------------------
## 2) 距离矩阵（CLR后欧氏距离= Aitchison距离）
## -----------------------------
dist_mat <- dist(clr, method = "euclidean")

## -----------------------------
## 3) 评估最佳k（CH + silhouette）
## -----------------------------
library(cluster)
library(tibble)

# 计算 Calinski–Harabasz 指数（在欧氏空间：用 CLR 后矩阵最合适）
calinski_harabasz <- function(x, cl){
  # x: n x p（样本×变量），cl: 长度 n 的簇编号
  x <- as.matrix(x)
  cl <- as.integer(cl)
  n <- nrow(x)
  k <- length(unique(cl))
  if (k < 2 || k >= n) return(NA_real_)
  
  grand <- colMeans(x)
  # 总离差
  tss <- sum(rowSums((x - matrix(grand, n, ncol(x), byrow = TRUE))^2))
  
  # 组内离差 WSS
  wss <- 0
  for (g in unique(cl)) {
    idx <- which(cl == g)
    cg <- colMeans(x[idx, , drop = FALSE])
    wss <- wss + sum(rowSums((x[idx, , drop = FALSE] - matrix(cg, length(idx), ncol(x), byrow = TRUE))^2))
  }
  
  bss <- tss - wss
  ch <- (bss / (k - 1)) / (wss / (n - k))
  ch
}

evaluate_k <- function(x_clr, dist_matrix, k_min = 2, k_max = 8){
  ks <- k_min:k_max
  
  ch <- sapply(ks, function(k){
    pam_fit <- pam(dist_matrix, k = k, diss = TRUE)
    calinski_harabasz(x_clr, pam_fit$clustering)
  })
  
  sil <- sapply(ks, function(k){
    pam_fit <- pam(dist_matrix, k = k, diss = TRUE)
    mean(silhouette(pam_fit$clustering, dist_matrix)[, 3])
  })
  
  tibble(k = ks, CH = ch, Silhouette = sil)
}


eval_df <- evaluate_k(clr, dist_mat, k_min = 2, k_max = 8)


p_ch <- ggplot(eval_df, aes(k, CH)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  theme_minimal(base_size = 12) +
  labs(title = "Calinski–Harabasz Index", x = "k", y = "CH (higher=better)")

p_ch

p_sil <- ggplot(eval_df, aes(k, Silhouette)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  theme_minimal(base_size = 12) +
  labs(title = "Silhouette Score", x = "k", y = "Silhouette (higher=better)")

p_sil

ggsave("./Figures/VFA_cluster_k_CH.pdf", p_ch, width = 4.5, height = 4)
ggsave("./Figures/VFA_cluster_k_Silhouette.pdf", p_sil, width = 4.5, height = 4)

## 手动选择k（建议看评估图后定；这里默认优先选 Silhouette 最大的k）
k_best <- eval_df$k[which.max(eval_df$Silhouette)]
message("Auto-selected k_best = ", k_best)

## -----------------------------
## 4) PAM 聚类
## -----------------------------
pam_res <- pam(dist_mat, k = 2, diss = TRUE)
cluster_raw <- pam_res$clustering  # 1..k

table(cluster_raw)

## -----------------------------
## 5) 按A:P给簇命名：Propionate-type / Acetate-type / Mixed...
## -----------------------------
ap <- derived$AP
names(ap) <- derived$SampleID

cluster_median_ap <- tapply(ap[names(cluster_raw)], cluster_raw, median)

# 从低到高：低AP = 丙酸型，高AP = 乙酸型
ord <- order(cluster_median_ap)

cluster_label_map <- rep("Mixed", length(cluster_median_ap))
names(cluster_label_map) <- names(cluster_median_ap)

cluster_label_map[names(cluster_median_ap)[ord[1]]] <- "Propionate-type"
cluster_label_map[names(cluster_median_ap)[ord[length(ord)]]] <- "Acetate-type"
# 其余保持Mixed（若k=2则只有两类；若k>3则中间全部Mixed）

rumenotype <- cluster_label_map[as.character(cluster_raw)]
rumenotype <- factor(rumenotype, levels = c("Propionate-type","Mixed","Acetate-type"))

table(rumenotype)

## -----------------------------
## 7) 导出结果
## -----------------------------
assignments <- data.frame(
  SampleID = rownames(clr),
  k = k_best,
  ClusterID = paste0("C", cluster_raw),
  Rumenotype = as.character(rumenotype),
  stringsAsFactors = FALSE
) %>%
  left_join(derived, by = "SampleID")

write.xlsx(assignments, "VFA_rumenotype_assignments.xlsx", overwrite = TRUE)

# 也把VFA比例、CLR矩阵、派生指标一起存下来，后面接宏基因组/非靶更方便
write.xlsx(list(
  "VFA_concentration" = data.frame(SampleID = raw[[id_col]], raw[, vfa_cols]),
  "VFA_proportion"    = data.frame(SampleID = rownames(prop), prop),
  "VFA_CLR"           = data.frame(SampleID = rownames(clr), clr),
  "Derived"           = derived,
  "k_evaluation"      = eval_df
), "VFA_features_and_k_eval.xlsx", overwrite = TRUE)

message("Done! 输出文件：",
        "\n- VFA_rumenotype_assignments.xlsx",
        "\n- VFA_features_and_k_eval.xlsx",
        "\n- VFA_cluster_k_CH.pdf",
        "\n- VFA_cluster_k_Silhouette.pdf",
        "\n- VFA_PCoA_rumenotypes.pdf")





# 取每个类型的样本ID
id1 <- assignments$SampleID[assignments$Rumenotype == "Acetate-type"]
id2 <- assignments$SampleID[assignments$Rumenotype == "Propionate-type"]

# 子矩阵（样本×VFA）
cluster1 <- clr[id1, , drop = FALSE]
cluster2 <- clr[id2, , drop = FALSE]

# 欧氏距离（CLR空间）
d1 <- dist(cluster1, method = "euclidean")
d2 <- dist(cluster2, method = "euclidean")

hc_1 <- hclust(d1, method = "average")
hc_2 <- hclust(d2, method = "average")

# 树图顺序：从左到右的 SampleID
ord1_id <- hc_1$labels[hc_1$order]
ord2_id <- hc_2$labels[hc_2$order]

# ✅ 用样本名重排 prop2（而不是用整数order）
prop2_ord <- prop2[c(ord1_id, ord2_id),]

# 画树（可选）
plot(hc_1, sub = "", xlab = "", hang = -1)
plot(hc_2, sub = "", xlab = "", hang = -1)

# 检查是否对齐
stopifnot(identical(rownames(prop2_ord)[1:length(ord1_id)], ord1_id))
stopifnot(identical(rownames(prop2_ord)[(length(ord1_id)+1):nrow(prop2_ord)], ord2_id))



prop2_ord_long <- prop2_ord %>% as.data.frame() %>% rownames_to_column(var = 'SampleID')

prop2_ord_long <- pivot_longer(prop2_ord_long,cols = -1,names_to = 'VFA',values_to = 'Propotion')


prop2_ord_long$SampleID <- factor(prop2_ord_long$SampleID, levels = unique(prop2_ord_long$SampleID))

prop2_ord_long$VFA <- factor(prop2_ord_long$VFA, levels = c("Acetic_acid","Propionic_acid","Butyric_acid","Isobutyric_acid",
                                                            "Valeric_acid","Isovaleric_acid"))

ggplot(data = prop2_ord_long, aes(SampleID,Propotion,fill = VFA))+  
  geom_bar(stat = "identity", position = "fill") +    
  labs(y = "Propotion",  
       fill = "VFA") + 
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = color)

ggsave("./Figures/VFA relative abundance.pdf", width = 10, height = 5)


# CQseq <- rbind(genus_CQseq[selected_genus,],apply(genus_CQseq_others,2,sum))
# WZ <- rbind(genus_WZ[selected_genus,],apply(genus_WZ_others,2,sum))
# 
# combined <- cbind(CQ,CQseq,WZ)

# 列注释

# annotation_col <- data.frame(row.names = c(metadata_CQ$Sample_ID,metadata_CQseq$SampleID,metadata_WZ$SampleID),
#                              OMT=c(metadata_CQ$OMT,metadata_CQseq$OMT,metadata_WZ$OMT),
#                              Group=c(metadata_CQ$Group,metadata_CQseq$Group,metadata_WZ$Group),
#                              Cohort=c(rep('CQ',nrow(metadata_CQ)),rep('CQseq',nrow(metadata_CQseq)),rep('WZ',nrow(metadata_WZ))),
#                              Fasting=c(metadata_CQ$Fasting_Glucose,rep(NA,length(metadata_CQseq$SampleID)),metadata_WZ$Fasting),
#                              OGTT1h=c(metadata_CQ$OGTT1h_Glucose,rep(NA,length(metadata_CQseq$SampleID)),metadata_WZ$GTT2h),
#                              OGTT2h=c(metadata_CQ$OGTT2h_Glucose,rep(NA,length(metadata_CQseq$SampleID)),metadata_WZ$GTT2h)
#                              )

annotation_col <- metadata_CQ %>% column_to_rownames(var = 'Sample_ID')

annotation_col['C_1_G0835','BMI'] <- NA

annotation <- annotation_col[,c(14,1:2,4,7,8:10)]

annotation$Fasting_Glucose <- ifelse(annotation$Fasting_Glucose>=5.1,"Yes",'NO')
annotation$OGTT1h_Glucose <- ifelse(annotation$OGTT1h_Glucose>=10,"Yes",'NO')
annotation$OGTT2h_Glucose <- ifelse(annotation$OGTT2h_Glucose>=8.5,"Yes",'NO')

pheatmap::pheatmap(CQ,annotation_col = annotation,cluster_cols = F, cluster_rows = F)

# 对AP比例进行统计

ggplot(assignments, aes(x=Rumenotype,y=AP, fill = Rumenotype, color = Rumenotype))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=2, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('Acetate-type'=darken('#2F5597',0.5),
                                'Propionate-type'=darken('#C00000',0.5)))+
  scale_fill_manual(values = c('Acetate-type'='#2F5597',
                               'Propionate-type'='#C00000'))+
  geom_signif(comparisons = list(c('Acetate-type','Propionate-type')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/AP.pdf',height = 4, width = 3)

assignments <- read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx",
                         sheet = 3)

library(ggplot2)
library(GGally)

vars <- c("Acetic_acid",
          "Propionic_acid",
          "Isobutyric_acid",
          "Butyric_acid",
          "Isovaleric_acid",
          "Valeric_acid")

p <- ggpairs(
  assignments,
  columns = vars,
  
  upper = list(
    continuous = wrap(
      "cor",
      method = "pearson",
      size = 4,
      color = "#C00000"   # 相关系数颜色
    )
  ),
  
  lower = list(
    continuous = wrap(
      "smooth",
      method = "lm",
      se = TRUE,
      color = "#C00000",     # 回归线
      fill = "#F4A6A6",      # CI区域
      alpha = 0.4,
      size = 0.4,
      point.colour = "#2F5597",  # 散点颜色
      point.size = 1.5
    )
  ),
  
  diag = list(
    continuous = wrap(
      "densityDiag",
      fill = "#4DBBD5",
      alpha = 0.7,
      color = "#2F5597"
    )
  )
) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # 去掉主网格
    panel.grid.minor = element_blank(),   # 去掉次网格
    strip.background = element_rect(fill="white"),
    strip.text = element_text(size = 11, face = "bold")
  )

p

ggsave(
  "./Figures/VFA_pairwise_regression_matrix_color.pdf",
  p,
  height = 8,
  width = 8
)

ggplot(assignments, aes(x=Rumenotype,y=Butyric_acid, fill = Rumenotype, color = Rumenotype))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=2, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        # axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('AEE'=darken('#2F5597',0.5),
                                'PEE'=darken('#C00000',0.5)))+
  scale_fill_manual(values = c('AEE'='#2F5597',
                               'PEE'='#C00000'))+
  geom_signif(comparisons = list(c('AEE','PEE')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/Butyric_acid.pdf',height = 4, width = 3)

library(ggplot2)
library(ggpubr)

p <- ggplot(assignments, aes(x = Acetic_acid, y = Butyric_acid)) +
  
  geom_point(
    size = 2.5,
    color = "#2C7FB8",
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    color = "#D95F02",
    fill = "#FDB863",
    se = TRUE,
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x = min(assignments$Acetic_acid),
    label.y = max(assignments$Butyric_acid),
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 4
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks.length = unit(2,"mm"),
    axis.title = element_text(size = 12),
    # axis.text = element_blank()
  ) +
  
  labs(
    x = "Acetic acid",
    y = "Butyric acid"
  )

p

ggsave(
  "./Figures/Butyric_acid_regression.pdf",
  p,
  height = 4,
  width = 4
)


p_Propionic_acid <- ggplot(assignments, aes(x = Propionic_acid, y = Butyric_acid)) +
  
  geom_point(
    size = 2.5,
    color = "#2C7FB8",
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    color = "#D95F02",
    fill = "#FDB863",
    se = TRUE,
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x = min(assignments$Acetic_acid),
    label.y = max(assignments$Butyric_acid),
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 4
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks.length = unit(2,"mm"),
    axis.title = element_text(size = 12),
    # axis.text = element_blank()
  ) +
  
  labs(
    x = "Propionic acid",
    y = "Butyric acid"
  )

p_Propionic_acid

ggsave(
  "./Figures/Propionic_acid_regression.pdf",
  p_Propionic_acid,
  height = 4,
  width = 4
)



ggplot(assignments, aes(x=Rumenotype,y=BCVFA, fill = Rumenotype, color = Rumenotype))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=2, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('Acetate-type'=darken('#2F5597',0.5),
                                'Propionate-type'=darken('#C00000',0.5)))+
  scale_fill_manual(values = c('Acetate-type'='#2F5597',
                               'Propionate-type'='#C00000'))+
  geom_signif(comparisons = list(c('Acetate-type','Propionate-type')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/BCVFA.pdf',height = 4, width = 3)


ggplot(assignments, aes(x=Rumenotype,y=TVFA, fill = Rumenotype, color = Rumenotype))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=2, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('Acetate-type'=darken('#2F5597',0.5),
                                'Propionate-type'=darken('#C00000',0.5)))+
  scale_fill_manual(values = c('Acetate-type'='#2F5597',
                               'Propionate-type'='#C00000'))+
  geom_signif(comparisons = list(c('Acetate-type','Propionate-type')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/BCVFA.pdf',height = 4, width = 3)


######### 氨氮含量统计#####################

NH3N <- read.xlsx('./代谢组数据/氨态氮总表.xlsx')

assignments <- read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx",
                         sheet = 1)

NH3N <- left_join(NH3N, assignments, by='SampleID')

NH3N <- NH3N %>% na.exclude()

cor_mat <- cor(NH3N[,c(2,6:9)], method = "pearson")

corrplot::corrplot(cor_mat)



p_NH3N <- ggplot(NH3N, aes(x = Propionic_acid , y = NH3N)) +
  
  geom_point(
    size = 2.5,
    color = "#2C7FB8",
    alpha = 0.8
  ) +
  
  geom_smooth(
    method = "lm",
    color = "#D95F02",
    fill = "#FDB863",
    se = TRUE,
    linewidth = 1
  ) +
  
  stat_cor(
    method = "pearson",
    label.x = min(NH3N$Propionic_acid),
    label.y = max(NH3N$NH3N),
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 4
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks.length = unit(2,"mm"),
    axis.title = element_text(size = 12),
    # axis.text = element_blank()
  ) +
  
  labs(
    x = "AP",
    y = "NH3N"
  )

p_NH3N

ggsave(
  "./Figures/Propionic_acid_regression.pdf",
  p_Propionic_acid,
  height = 4,
  width = 4
)
  
ggplot(NH3N, aes(x=Rumenotype,y=NH3N, fill = Rumenotype, color = Rumenotype))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=2, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        # axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('Acetate-type'=darken('#2F5597',0.5),
                                'Propionate-type'=darken('#C00000',0.5)))+
  scale_fill_manual(values = c('Acetate-type'='#2F5597',
                               'Propionate-type'='#C00000'))+
  geom_signif(comparisons = list(c('Acetate-type','Propionate-type')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/NH3N.pdf',height = 4, width = 3)



# 使用PHATE降维分析

# 安装并加载reticulate包
library(reticulate)
library(phateR)

Sys.setenv(RETICULATE_PYTHON = "C:/Users/sheng/anaconda3/envs/r-phate/python.exe")
use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)

phate <- import("phate")

# 使用预先计算的距离矩阵，设置knn.dist.method为"precomputed"
phate_result <- phate(
  dist_mat,
  knn.dist.method = "precomputed",  # 指定使用预计算距离
  knn = 5,                          # 近邻数，根据数据调整
  ndim = 3,                         # 输出维度，通常2或3
  t = "auto",                       # 自动选择扩散时间
  n.landmark = 2000,                # 使用地标点加速计算，可选
  gamma = 1,                        # 熵正则化参数
  verbose = TRUE                    # 显示进度信息
)

# 获取PHATE坐标
phate_coords <- phate_result$embedding %>% as.data.frame() %>% rownames_to_column(var = 'SampleID')

phate_coords <- left_join(phate_coords, assignments,by='SampleID')

colors <- c('Acetate-type'='#2F5597',
            'Propionate-type'='#C00000')


# 可视化
ggplot(phate_coords, aes(PHATE1, PHATE2, color = Rumenotype)) +
  geom_point(size=2) +
  theme_bw()+
  theme(axis.title = element_text(size = 10),
        axis.title.y = element_text(angle = 90),
        axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        legend.position = 'inside',
        legend.position.inside = c(0.5,0.8),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = colors)

ggsave('./Figures/代谢型PHATE降维-2d.pdf',width = 4, height = 4)


ph <- assignments[, c("SampleID","AP")]
ph$FID <- ph$SampleID
ph$IID <- ph$SampleID
ph <- ph[, c("FID","IID","AP")]
write.table(ph, "./Figures/AP_pheno.tsv", sep="\t", row.names=FALSE, quote=FALSE)


################################
# Bootstrap一致性代码（推荐 subsampling 80% × 200次）
################################
suppressPackageStartupMessages({
  library(cluster)        # pam
  library(mclust)         # adjustedRandIndex
  library(aricode)        # NMI (optional)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

## -----------------------------
## 0) 输入：你的全量结果
## -----------------------------
# clr: n x p, rownames=SampleID
# cluster_raw: 全量PAM聚类结果（长度n，names=SampleID）
stopifnot(!is.null(rownames(clr)))
stopifnot(!is.null(names(cluster_raw)))
stopifnot(all(rownames(clr) %in% names(cluster_raw)))

# 全量标签（按 clr 行顺序对齐）
full_ids <- rownames(clr)
full_lab <- as.integer(cluster_raw[full_ids])

## -----------------------------
## 1) 一个小工具：把“子集聚类标签”对齐到“全量标签”的簇编号
##    （解决 PAM 的簇编号在每次重跑可能互换的问题）
## -----------------------------
align_labels_to_full <- function(full_labels_sub, sub_labels){
  # full_labels_sub: 子集样本在全量聚类里的标签（长度 m）
  # sub_labels:      子集重新聚类得到的标签（长度 m）
  # 输出：映射后的 sub_labels，使其尽可能与 full_labels_sub 对齐
  tab <- table(sub_labels, full_labels_sub)
  
  # 贪心匹配：每个sub簇映射到最多重叠的full簇
  map <- apply(tab, 1, function(x) as.integer(names(which.max(x))))
  aligned <- map[as.character(sub_labels)]
  as.integer(aligned)
}

## -----------------------------
## 2) Bootstrap / Subsampling 稳定性评估
## -----------------------------
bootstrap_cluster_stability <- function(clr, full_lab, frac = 0.8, B = 200, k = length(unique(full_lab)),
                                        seed = 1){
  set.seed(seed)
  ids <- rownames(clr)
  n <- length(ids)
  m <- max(2, floor(n * frac))
  
  # 记录每次一致性
  metrics <- vector("list", B)
  
  # 记录每个样本：被抽到次数 + 与全量一致次数
  hit <- setNames(integer(n), ids)
  agree <- setNames(integer(n), ids)
  
  for(b in seq_len(B)){
    sub_ids <- sample(ids, size = m, replace = FALSE)
    sub_x <- clr[sub_ids, , drop = FALSE]
    
    # 子集聚类
    sub_dist <- dist(sub_x, method = "euclidean")
    pam_fit <- pam(sub_dist, k = k, diss = TRUE)
    sub_lab <- pam_fit$clustering
    
    # 子集样本在全量中的标签
    full_lab_sub <- full_lab[match(sub_ids, ids)]
    
    # 对齐簇编号（避免1/2互换导致假不一致）
    sub_lab_aligned <- align_labels_to_full(full_lab_sub, sub_lab)
    
    # 指标：与全量标签的一致性（子集内）
    ari <- mclust::adjustedRandIndex(full_lab_sub, sub_lab_aligned)
    
    # NMI（可选，若没装aricode就跳过）
    nmi <- NA_real_
    if (requireNamespace("aricode", quietly = TRUE)) {
      nmi <- aricode::NMI(full_lab_sub, sub_lab_aligned)
    }
    
    metrics[[b]] <- tibble(iter = b, ARI = ari, NMI = nmi, n_sub = m)
    
    # 更新样本级稳定性
    hit[sub_ids] <- hit[sub_ids] + 1L
    agree[sub_ids] <- agree[sub_ids] + as.integer(sub_lab_aligned == full_lab_sub)
  }
  
  metrics_df <- bind_rows(metrics)
  
  sample_stab <- tibble(
    SampleID = ids,
    n_hit = as.integer(hit),
    n_agree = as.integer(agree),
    stability = ifelse(hit > 0, agree / hit, NA_real_)
  )
  
  list(metrics = metrics_df, sample_stability = sample_stab)
}

res_stab <- bootstrap_cluster_stability(
  clr = clr,
  full_lab = full_lab,
  frac = 0.8,
  B = 200,
  k = length(unique(full_lab)),
  seed = 123
)

metrics_df <- res_stab$metrics
sample_stab <- res_stab$sample_stability

## -----------------------------
## 3) 可视化：ARI/NMI分布 + 样本稳定性
## -----------------------------
p_ari <- ggplot(metrics_df, aes(x = ARI)) +
  geom_histogram(bins = 30) +
  theme_minimal(base_size = 12) +
  labs(title = "Bootstrap stability (subsampling 80%)", x = "ARI vs full-data clustering", y = "Count")

p_nmi <- ggplot(metrics_df, aes(x = NMI)) +
  geom_histogram(bins = 30) +
  theme_minimal(base_size = 12) +
  labs(title = "Bootstrap stability (subsampling 80%)", x = "NMI vs full-data clustering", y = "Count")

p_sample <- ggplot(sample_stab, aes(x = stability)) +
  geom_histogram(bins = 30) +
  theme_minimal(base_size = 12) +
  labs(title = "Per-sample assignment stability", x = "Agreement probability", y = "Count")

print(p_ari)
if (!all(is.na(metrics_df$NMI))) print(p_nmi)
print(p_sample)

## 导出结果
dir.create("./Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("./Figures/VFA_cluster_bootstrap_ARI.pdf", p_ari, width = 4.5, height = 4)
if (!all(is.na(metrics_df$NMI))) ggsave("./Figures/VFA_cluster_bootstrap_NMI.pdf", p_nmi, width = 4.5, height = 4)
ggsave("./Figures/VFA_cluster_sample_stability.pdf", p_sample, width = 4.5, height = 4)

write.table(metrics_df, "./Figures/VFA_cluster_bootstrap_metrics.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(sample_stab, "./Figures/VFA_cluster_sample_stability.tsv", sep="\t", quote=FALSE, row.names=FALSE)

summary(metrics_df$ARI)


###########################
# k=2 vs k=3 的 bootstrap 稳定性对比图”（ARI/NMI）
###########################

suppressPackageStartupMessages({
  library(cluster)
  library(mclust)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# 对齐簇编号（避免标签互换）
align_labels_to_ref <- function(ref, lab){
  tab <- table(lab, ref)
  map <- apply(tab, 1, function(x) as.integer(names(which.max(x))))
  as.integer(map[as.character(lab)])
}

bootstrap_stability_for_k <- function(clr, k, frac = 0.8, B = 200, seed = 1){
  set.seed(seed)
  ids <- rownames(clr)
  n <- length(ids)
  m <- max(2, floor(n * frac))
  
  # 全量聚类作为 reference
  d_full <- dist(clr, method = "euclidean")
  ref <- pam(d_full, k = k, diss = TRUE)$clustering
  names(ref) <- ids
  
  ari_vec <- numeric(B)
  
  for (b in seq_len(B)){
    sub_ids <- sample(ids, m, replace = FALSE)
    sub_x <- clr[sub_ids, , drop = FALSE]
    d_sub <- dist(sub_x, method = "euclidean")
    lab_sub <- pam(d_sub, k = k, diss = TRUE)$clustering
    
    ref_sub <- ref[sub_ids]
    lab_sub2 <- align_labels_to_ref(ref_sub, lab_sub)
    
    ari_vec[b] <- mclust::adjustedRandIndex(ref_sub, lab_sub2)
  }
  
  tibble(iter = 1:B, k = k, ARI = ari_vec)
}

res_k2 <- bootstrap_stability_for_k(clr, k = 2, frac = 0.8, B = 200, seed = 123)
res_k3 <- bootstrap_stability_for_k(clr, k = 3, frac = 0.8, B = 200, seed = 123)
res_all <- bind_rows(res_k2, res_k3)

# 打印摘要
print(res_all %>% group_by(k) %>%
        summarise(n = n(), median_ARI = median(ARI), mean_ARI = mean(ARI),
                  q1 = quantile(ARI, 0.25), q3 = quantile(ARI, 0.75)))

# 画对比图
p <- ggplot(res_all, aes(x = ARI)) +
  geom_histogram(bins = 30) +
  facet_wrap(~k, nrow = 1, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Subsampling stability comparison", x = "ARI vs full-data clustering", y = "Count")

dir.create("./Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("./Figures/VFA_k2_vs_k3_bootstrap_ARI.pdf", p, width = 8.5, height = 4)
p

#######################
#
#######################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cluster)
  library(scales)
})

## 1) 计算 k=2 / k=3 聚类（如果你已有结果可跳过）
dist_mat <- dist(clr, method = "euclidean")

pam2 <- pam(dist_mat, k = 2, diss = TRUE)
pam3 <- pam(dist_mat, k = 3, diss = TRUE)

cluster2 <- pam2$clustering
cluster3 <- pam3$clustering

## 2) 汇总簇大小
tab2 <- as.data.frame(table(cluster2), stringsAsFactors = FALSE) %>%
  transmute(k = "k=2", Cluster = paste0("C", cluster2), n = Freq)

tab3 <- as.data.frame(table(cluster3), stringsAsFactors = FALSE) %>%
  transmute(k = "k=3", Cluster = paste0("C", cluster3), n = Freq)

size_df <- bind_rows(tab2, tab3) %>%
  group_by(k) %>%
  mutate(
    frac = n / sum(n),
    lab = paste0(Cluster, "\n", "n=", n, "\n(", percent(frac, accuracy = 0.1), ")"),
    # ✅ 突出：k=3 中最小的簇
    is_tiny = (k == "k=3" & n == min(n))
  ) %>%
  ungroup()

## 3) 画 donut（用极坐标 + 空心）
p_donut <- ggplot(size_df, aes(x = 2, y = frac, fill = Cluster)) +
  geom_col(color = "white", linewidth = 0.4) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +   # 控制 donut 厚度
  facet_wrap(~k, nrow = 1) +
  theme_void(base_size = 12) +
  guides(fill = guide_legend(title = "Cluster")) +
  # ✅ 用文字标出小簇（只给 tiny 那块打标签）
  geom_text(
    data = size_df %>% filter(is_tiny),
    aes(label = paste0("Tiny cluster\n", "n=", n, "\n", percent(frac, accuracy = 0.01))),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  labs(title = "Cluster size distribution (k=2 vs k=3)")

dir.create("./Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("./Figures/VFA_cluster_size_donut_k2_k3.pdf", p_donut, width = 9, height = 4.2)

p_donut

####
#VFA 平均相对丰度
####

apply(prop2, 2, mean)


vfa_prop <- read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 1-VFA_rumenotype_assignments.xlsx",
                      sheet = 3)


library(ggtern)

ggplot(vfa_prop, aes(x = Rumenotype, y =  Acetic_acid/Propionic_acid, fill = Rumenotype)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("AEE"="#2F5597","PEE"="#C00000")) +
  theme_bw()



