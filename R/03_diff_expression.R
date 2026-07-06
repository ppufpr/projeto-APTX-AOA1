# ============================================================
# 03_diff_expression.R — Expressão diferencial (DESeq2 fatorial)
# Contrastes de interesse:
#  (a) EFEITO DO KO        : APTX-KO vs WT (genótipo)        -> baseline/referência
#  (b) EFEITO DO ESTÍMULO  : IS vs NS
#  (c) INTERAÇÃO           : o KO muda a resposta ao estímulo?
# ============================================================
source("R/config.R")
suppressMessages({ library(DESeq2); library(data.table) })

dds <- readRDS(file.path(DIR_PROCESSED, "dds.rds"))
message("Coeficientes disponíveis: ", paste(resultsNames(dds), collapse = ", "))

save_res <- function(res, tag) {
  df <- as.data.frame(res); df$gene <- rownames(df)
  df <- df[order(df$padj), ]
  fwrite(df, file.path(DIR_TABLES, paste0("de_", tag, "_all.csv")))
  degs <- subset(df, !is.na(padj) & padj < DE_PADJ & abs(log2FoldChange) >= DE_LFC)
  fwrite(degs, file.path(DIR_TABLES, paste0("de_", tag, "_degs.csv")))
  message("DE [", tag, "]: ", nrow(degs), " DEGs (padj<", DE_PADJ, ", |LFC|>=", DE_LFC, ")")
  invisible(df)
}

# (a) Efeito do genótipo (KO vs WT) — contraste principal
res_geno <- results(dds, contrast = c("genotype", "KO", "WT"))
save_res(res_geno, "genotype_KO_vs_WT")

# (b) Efeito do estímulo (IS vs NS)
res_stim <- results(dds, contrast = c("stimulation", "IS", "NS"))
save_res(res_stim, "stimulation_IS_vs_NS")

# (c) Interação genótipo x estímulo (último termo do desenho)
int_name <- tail(resultsNames(dds), 1)
res_int  <- results(dds, name = int_name)
save_res(res_int, "interaction")

# Volcano do contraste principal (genótipo)
df <- as.data.frame(res_geno)
png(file.path(DIR_FIGURES, "volcano_genotype.png"), width = 1200, height = 900, res = 150)
with(df, plot(log2FoldChange, -log10(pmax(padj, 1e-300)),
              pch = 20, col = ifelse(!is.na(padj) & padj < DE_PADJ, "firebrick", "grey60"),
              xlab = "log2 FC (APTX-KO / WT)", ylab = "-log10(FDR)",
              main = "Expressão diferencial — APTX-KO x WT"))
abline(h = -log10(DE_PADJ), lty = 2); abline(v = c(-DE_LFC, DE_LFC), lty = 2)
dev.off()

message("Expressão diferencial concluída (3 contrastes salvos em results/tables/).")
