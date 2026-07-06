# ============================================================
# 01_download.R — baixa os dados do GEO (GSE245766)
# ============================================================
source("R/config.R")
suppressMessages({ library(GEOquery); library(data.table) })

# ---- 1. Arquivos suplementares (matriz de contagens) ----
message("Baixando arquivos suplementares do ", GSE, " ...")
getGEOSuppFiles(GSE, baseDir = DIR_RAW, makeDirectory = TRUE)

# ---- 2. Metadados (phenoData) para mapear amostra -> condição ----
message("Baixando metadados da série ...")
gse_meta <- getGEO(GSE, GSEMatrix = TRUE, getGPL = FALSE)
pheno <- as.data.frame(Biobase::pData(gse_meta[[1]]))
fwrite(pheno, file.path(DIR_PROCESSED, "pheno.csv"))

message("\nArquivos baixados:")
print(list.files(DIR_RAW, recursive = TRUE))

message("\n>>> TODO: abra data/processed/pheno.csv e identifique as colunas que indicam ",
        "GENÓTIPO (APTX-KO vs WT) e ESTÍMULO (IS vs NS). Use-as em 02_preprocess.R.")
