


tsne <- read_tsv("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/微生物/MAG_COG_cluster/Tables/MAG_embeddings_clusters.tsv")

diffmag <- readxl::read_excel("D:/BaiduSyncdisk/内师大/在研课题/2025宿主与微生物定植/FinalFigure/Supplemental file 4-AEE PEE differential MAGs.xlsx", sheet = 1) %>% as.data.frame()

diffmag <- left_join(diffmag, tsne, by= 'MAG_formal')

library(vegan)

Y <- as.matrix(diffmag[, c("TSNE1", "TSNE2")])

adon_res <- adonis2(
  Y ~ Rumenotype,
  data = diffmag,
  permutations = 999,
  method = "euclidean"
)

print(adon_res)
summary(adon_res)
