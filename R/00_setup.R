# ============================================================
# 00_setup.R — instala e carrega as dependências
# Rodar uma única vez antes do pipeline.
# ============================================================

## --- Pacotes CRAN ---
cran_pkgs <- c(
  "data.table",   # leitura rápida de tabelas
  "matrixStats",  # estatísticas por linha/coluna
  "igraph",       # redes / centralidade
  "WGCNA",        # redes de co-expressão
  "CoDiNA",       # co-expressão diferencial (pacote da Profa. Deisy Gysi)
  "glmnet",       # classificador regularizado
  "caret",        # validação cruzada
  "pROC",         # curvas ROC / AUC
  "ggplot2",      # gráficos
  "pheatmap"      # heatmaps
)

## --- Pacotes Bioconductor ---
bioc_pkgs <- c(
  "GEOquery",        # download do GEO
  "DESeq2",          # normalização + expressão diferencial
  "limma",           # apoio
  "edgeR",           # CPM / filtro
  "clusterProfiler", # enriquecimento GO/KEGG
  "org.Hs.eg.db",    # anotação HUMANO
  "AnnotationDbi",
  "STRINGdb"         # rede de interações (Medicina de Redes)
)

## Instalador CRAN
new_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(new_cran)) install.packages(new_cran, repos = "https://cloud.r-project.org")

## Instalador Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
new_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

## Carrega tudo
invisible(lapply(c(cran_pkgs, bioc_pkgs), function(p)
  suppressMessages(library(p, character.only = TRUE))))

message("Setup concluído. Pacotes carregados: ",
        paste(c(cran_pkgs, bioc_pkgs), collapse = ", "))
