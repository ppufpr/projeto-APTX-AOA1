# ============================================================
# 06_enrichment.R — Enriquecimento funcional GO/KEGG (humano)
# Interpreta os DEGs do efeito do KO (e/ou módulos WGCNA).
# Espera-se destaque para: reparo de DNA, resposta imune inata
# (cGAS-STING, RIG-I/MAVS), sinalização antiviral.
# ============================================================
source("R/config.R")
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi) })

deg_file <- file.path(DIR_TABLES, "de_genotype_KO_vs_WT_degs.csv")
stopifnot(file.exists(deg_file))
genes <- read.csv(deg_file)$gene
genes <- sub("\\.\\d+$", "", genes)   # remove sufixo de versão Ensembl (ENSG...->.N)

# Mapear ID -> ENTREZ (ajuste keytype ao ID do dataset: SYMBOL ou ENSEMBL).
# tryCatch: mapIds lança erro (não NA) quando NENHUMA chave é válida -> permite fallback.
id_map <- function(keytype)
  tryCatch(
    suppressWarnings(AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes, keytype = keytype,
                                           column = "ENTREZID", multiVals = "first")),
    error = function(e) setNames(rep(NA_character_, length(genes)), genes))
entrez <- id_map("SYMBOL")
if (mean(is.na(entrez)) > 0.8) entrez <- id_map("ENSEMBL")
entrez <- unique(na.omit(entrez))
message("Genes mapeados para ENTREZ: ", length(entrez))

# GO (Processo Biológico)
ego <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ, readable = TRUE)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write.csv(as.data.frame(ego), file.path(DIR_TABLES, "enrich_GO_BP.csv"), row.names = FALSE)
  png(file.path(DIR_FIGURES, "enrich_GO_BP.png"), width = 1400, height = 1000, res = 140)
  print(dotplot(ego, showCategory = 20) + ggplot2::ggtitle("GO BP — APTX-KO x WT"))
  dev.off()
}

# KEGG
ekegg <- enrichKEGG(entrez, organism = KEGG_ORG, pAdjustMethod = "BH", pvalueCutoff = ENRICH_PADJ)
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  write.csv(as.data.frame(ekegg), file.path(DIR_TABLES, "enrich_KEGG.csv"), row.names = FALSE)
  png(file.path(DIR_FIGURES, "enrich_KEGG.png"), width = 1400, height = 1000, res = 140)
  print(dotplot(ekegg, showCategory = 20) + ggplot2::ggtitle("KEGG — APTX-KO x WT"))
  dev.off()
}

message("Enriquecimento concluído. Procure por: DNA repair, innate immune response, ",
        "cytosolic DNA-sensing (cGAS-STING), RIG-I-like receptor signaling.")
