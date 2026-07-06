# ============================================================
# 04_wgcna.R — Rede de co-expressão e detecção de módulos (WGCNA)
# Obs.: n = 12 (2x2). WGCNA aqui é exploratório — interprete os módulos
# com cautela; foco em padrões robustos, não em arestas isoladas.
# ============================================================
source("R/config.R")
suppressMessages({ library(WGCNA); library(matrixStats) })
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

vsd     <- readRDS(file.path(DIR_PROCESSED, "vsd.rds"))
coldata <- readRDS(file.path(DIR_PROCESSED, "coldata.rds"))

# WGCNA: amostras nas LINHAS, genes nas COLUNAS. Top genes mais variáveis.
n_top <- min(5000, nrow(vsd))
top_genes <- order(rowVars(vsd), decreasing = TRUE)[seq_len(n_top)]
datExpr <- t(vsd[top_genes, ])

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# Soft-threshold (beta)
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)

# Power recomendado pelo WGCNA quando o ajuste scale-free falha (comum com n pequeno).
# Tabela oficial p/ rede SIGNED: <20 amostras -> 18; 20-30 -> 16; 30-40 -> 14; >40 -> 12.
recommended_power <- function(n_samples) {
  if (n_samples < 20) 18L else if (n_samples < 30) 16L else if (n_samples < 40) 14L else 12L
}
# Aceita a estimativa do pickSoftThreshold só se for plausível p/ signed (>=12);
# caso contrário (NA ou beta baixo espúrio, como o beta=3 com n=12), usa a recomendação.
power <- if (!is.null(WGCNA_POWER)) {
  WGCNA_POWER
} else if (!is.na(sft$powerEstimate) && sft$powerEstimate >= 12L) {
  sft$powerEstimate
} else {
  recommended_power(nrow(datExpr))
}
message("Soft-threshold: estimado pelo fit = ",
        ifelse(is.na(sft$powerEstimate), "NA", sft$powerEstimate),
        " | beta usado = ", power,
        if (is.na(sft$powerEstimate) || sft$powerEstimate < 12L)
          " (fallback recomendado WGCNA: signed, n pequeno)" else "")

png(file.path(DIR_FIGURES, "wgcna_scalefree.png"), width = 1200, height = 600, res = 130)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="power", ylab="Scale-free R^2", type="b", main="Ajuste scale-free"); abline(h=0.8, col="red", lty=2)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="power", ylab="Conectividade média", type="b", main="Conectividade")
dev.off()

# Módulos
net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = WGCNA_MINMODSZ, mergeCutHeight = 0.25,
                        numericLabels = TRUE, saveTOMs = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
message("Módulos detectados: ", length(unique(moduleColors)))

# Associação módulo x traços (genótipo e estímulo)
traits <- data.frame(
  genotype = as.numeric(coldata[rownames(datExpr), "genotype"]),   # WT=1, KO=2
  stim     = as.numeric(coldata[rownames(datExpr), "stimulation"]) # NS=1, IS=2
)
mtc <- cor(MEs, traits, use = "p")
mtp <- corPvalueStudent(mtc, nrow(datExpr))
out <- data.frame(module = rownames(mtc),
                  cor_genotype = mtc[,"genotype"], p_genotype = mtp[,"genotype"],
                  cor_stim = mtc[,"stim"], p_stim = mtp[,"stim"])
out <- out[order(out$p_genotype), ]
write.csv(out, file.path(DIR_TABLES, "wgcna_module_trait.csv"), row.names = FALSE)

# Gene -> módulo
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
          file.path(DIR_TABLES, "wgcna_gene_module.csv"), row.names = FALSE)

saveRDS(list(net=net, moduleColors=moduleColors, MEs=MEs, datExpr=datExpr, power=power),
        file.path(DIR_PROCESSED, "wgcna.rds"))

# Em qual módulo cai o próprio APTX e os genes-semente?
gm <- setNames(moduleColors, colnames(datExpr))
seed_mod <- gm[intersect(SEED_GENES, names(gm))]
if (length(seed_mod)) {
  message("Módulo dos genes-semente:")
  print(seed_mod)
}

message("WGCNA concluído. Associação módulo-condição em results/tables/wgcna_module_trait.csv")
