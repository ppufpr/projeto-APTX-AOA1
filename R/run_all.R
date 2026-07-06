# ============================================================
# run_all.R — executa o pipeline na ordem correta
# Pré-requisito: ter rodado source("R/00_setup.R") uma vez.
# ============================================================
t0 <- Sys.time()

steps <- c(
  "R/01_download.R",          # baixa GSE245766
  "R/02_preprocess.R",        # metadados 2x2 + filtro + vst
  "R/03_diff_expression.R",   # DESeq2 fatorial (baseline/referência)
  "R/04_wgcna.R",             # rede de co-expressão + módulos
  "R/05_diff_coexpression.R", # CoDiNA (APTX-KO x WT)
  "R/06_enrichment.R",        # GO/KEGG (humano)
  "R/07_classification.R",    # métricas: AUC, precisão, recall, F1
  "R/08_network_medicine.R"   # módulo de reparo de DNA (STRING)
)

for (s in steps) {
  message("\n========== ", s, " ==========")
  source(s, echo = FALSE)
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("\nPipeline concluído em ",
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min.")
