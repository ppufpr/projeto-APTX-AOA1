# ============================================================
# 02_preprocess.R — contagens -> metadados (2x2), filtro e normalização
# ============================================================
source("R/config.R")
suppressMessages({ library(data.table); library(edgeR); library(DESeq2); library(matrixStats) })

# ------------------------------------------------------------
# 1. Carregar matriz de contagens (GSE245766)
#    O arquivo suplementar costuma ser um .txt/.csv.gz com genes nas linhas.
# ------------------------------------------------------------
cnt_path <- list.files(file.path(DIR_RAW, GSE), pattern = "count|matrix|raw",
                       full.names = TRUE, ignore.case = TRUE)
cnt_path <- cnt_path[grepl("\\.(txt|csv|tsv)(\\.gz)?$", cnt_path, ignore.case = TRUE)]
stopifnot(length(cnt_path) >= 1)
counts <- as.data.frame(fread(cnt_path[1]))
rownames(counts) <- counts[[1]]; counts[[1]] <- NULL
counts <- as.matrix(counts)
message("Contagens: ", nrow(counts), " genes x ", ncol(counts), " amostras")

# ------------------------------------------------------------
# 2. Definir condições (genótipo x estímulo)
#    TODO: ajuste as regras conforme data/processed/pheno.csv.
#    Espera-se 4 grupos: WT_NS, WT_IS, KO_NS, KO_IS (3 réplicas cada).
# ------------------------------------------------------------
sn <- colnames(counts)
genotype    <- factor(ifelse(grepl("KO|APTX|knock", sn, ignore.case = TRUE), "KO", "WT"),
                      levels = LEVELS_GENOTYPE)
# Estímulo: distinguir "ImmStim" (estimulado = IS) de "NoStim" (não estimulado = NS).
# CUIDADO: ambos os rótulos contêm "stim"; classificar primeiro o que é NS.
stimulation <- factor(ifelse(grepl("NoStim|unstim|NS\\b|_NS", sn, ignore.case = TRUE), "NS",
                      ifelse(grepl("ImmStim|stim|LPS|poly|IS\\b|_IS", sn, ignore.case = TRUE), "IS", "NS")),
                      levels = LEVELS_STIMULATION)
coldata <- data.frame(row.names = sn, genotype = genotype, stimulation = stimulation)
coldata$group <- factor(paste(coldata$genotype, coldata$stimulation, sep = "_"))
print(coldata)
stopifnot(nlevels(droplevels(genotype))    == 2,
          nlevels(droplevels(stimulation)) == 2)

# ------------------------------------------------------------
# 3. Filtro de baixa expressão
# ------------------------------------------------------------
dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = coldata$group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
message("Após filtro: ", nrow(dge), " genes")

# ------------------------------------------------------------
# 4. DESeq2 com desenho fatorial + normalização vst (para redes)
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(dge$counts),
                              colData   = coldata,
                              design    = ~ genotype + stimulation + genotype:stimulation)
dds <- DESeq(dds)
vsd <- assay(vst(dds, blind = FALSE))

# ------------------------------------------------------------
# 5. Salvar
# ------------------------------------------------------------
saveRDS(dds,     file.path(DIR_PROCESSED, "dds.rds"))
saveRDS(vsd,     file.path(DIR_PROCESSED, "vsd.rds"))
saveRDS(coldata, file.path(DIR_PROCESSED, "coldata.rds"))

message("Pré-processamento concluído (desenho 2x2).")
