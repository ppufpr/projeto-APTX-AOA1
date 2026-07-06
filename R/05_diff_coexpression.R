# ============================================================
# 05_diff_coexpression.R — Co-expressão diferencial (CoDiNA)
# Compara a rede de co-expressão APTX-KO vs WT e identifica ligações
# e genes cuja conectividade muda com a perda de APTX.
# ============================================================
source("R/config.R")
suppressMessages({ library(CoDiNA); library(matrixStats) })

vsd     <- readRDS(file.path(DIR_PROCESSED, "vsd.rds"))
coldata <- readRDS(file.path(DIR_PROCESSED, "coldata.rds"))

# Genes mais variáveis (rede tratável)
n_top <- min(2000, nrow(vsd))
top_genes <- order(rowVars(vsd), decreasing = TRUE)[seq_len(n_top)]
expr <- vsd[top_genes, ]

s_ko <- rownames(coldata)[coldata$genotype == "KO"]
s_wt <- rownames(coldata)[coldata$genotype == "WT"]

cor_ko <- cor(t(expr[, s_ko]), method = "spearman")
cor_wt <- cor(t(expr[, s_wt]), method = "spearman")

to_edges <- function(M, thr = 0.7) {
  M[lower.tri(M, diag = TRUE)] <- NA
  idx <- which(abs(M) >= thr, arr.ind = TRUE)
  data.frame(Node.1 = rownames(M)[idx[,1]], Node.2 = colnames(M)[idx[,2]], weight = M[idx])
}
edges_ko <- to_edges(cor_ko)
edges_wt <- to_edges(cor_wt)

dn <- MakeDiffNet(Data = list(edges_ko[,1:3], edges_wt[,1:3]), Code = c("KO", "WT"))
write.csv(as.data.frame(dn), file.path(DIR_TABLES, "codina_diffnet.csv"), row.names = FALSE)

nodes <- ClusterNodes(dn, cutoff.external = 0, cutoff.internal = 1)
write.csv(nodes, file.path(DIR_TABLES, "codina_nodes.csv"), row.names = FALSE)

# Onde está o APTX e os genes-semente na rede diferencial?
seed_in <- nodes[nodes$Node %in% SEED_GENES, , drop = FALSE]
if (nrow(seed_in)) {
  write.csv(seed_in, file.path(DIR_TABLES, "codina_seed_nodes.csv"), row.names = FALSE)
  message("Genes-semente na rede diferencial:"); print(seed_in)
}

message("Co-expressão diferencial (CoDiNA) concluída. Saídas em results/tables/codina_*.csv")
