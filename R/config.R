# ============================================================
# config.R — parâmetros centrais do projeto (APTX / AOA1)
# ============================================================

## Accession GEO
GSE <- "GSE245766"   # microglia humana HMC3; APTX-KO/WT x estímulo imune (2x2)

## Organismo (HUMANO)
ORG_DB     <- "org.Hs.eg.db"
KEGG_ORG   <- "hsa"
STRING_TAX <- 9606          # Homo sapiens (STRINGdb)

## Desenho experimental (fatorial 2x2)
#  genotype:    WT  | KO  (APTX knockout)
#  stimulation: NS  | IS  (not stimulated / immune stimulated)
LEVELS_GENOTYPE    <- c("WT", "KO")
LEVELS_STIMULATION <- c("NS", "IS")

## Genes-semente: módulo das ataxias de reparo de DNA (Medicina de Redes)
#  APTX=AOA1, SETX=AOA2, PNKP=AOA4, ATM=AT, TDP1=SCAN1; + parceiros de SSBR/NHEJ.
SEED_GENES <- c("APTX", "SETX", "PNKP", "ATM", "TDP1",
                "XRCC1", "XRCC4", "LIG4", "LIG3", "PARP1")

## Caminhos
DIR_RAW       <- "data/raw"
DIR_PROCESSED <- "data/processed"
DIR_TABLES    <- "results/tables"
DIR_FIGURES   <- "results/figures"

## Reprodutibilidade
SEED <- 42

## Parâmetros de análise
DE_PADJ        <- 0.05
DE_LFC         <- 1.0
WGCNA_POWER    <- NULL   # NULL => escolher pelo scale-free fit
WGCNA_MINMODSZ <- 30
ENRICH_PADJ    <- 0.05
STRING_SCORE   <- 400    # confiança mínima de aresta no STRING (0-1000)

## Cria diretórios de saída
for (d in c(DIR_RAW, DIR_PROCESSED, DIR_TABLES, DIR_FIGURES)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(SEED)
