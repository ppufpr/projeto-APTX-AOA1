# ============================================================
# 07_classification.R — Classificador APTX-KO x WT + MÉTRICAS
# Gera AUC, precisão, recall, F1 e matriz de confusão.
# Compara DOIS conjuntos de atributos:
#   (A) BASELINE: top DEGs (efeito do genótipo)  -> referência
#   (B) REDE:     eigengenes dos módulos WGCNA    -> proposto
# Validação cruzada leave-one-out (LOO) dado o n pequeno (n=12).
#
# RIGOR (evita viés de seleção; cf. Ambroise & McLachlan 2002; ESL cap. 7):
#   - Baseline (A): a SELEÇÃO de DEGs é refeita DENTRO de cada fold LOO
#     (limma na subamostra de treino, controlando o estímulo). Assim a
#     amostra de teste nunca influencia quais genes entram no modelo.
#   - Rede (B): os módulos/eigengenes do WGCNA são NÃO SUPERVISIONADOS
#     (co-expressão, sem olhar o genótipo) -> não há vazamento de rótulo;
#     usados como conjunto fixo de atributos.
# ============================================================
source("R/config.R")
suppressMessages({ library(glmnet); library(caret); library(pROC)
                   library(limma); library(matrixStats) })
set.seed(SEED)

vsd     <- readRDS(file.path(DIR_PROCESSED, "vsd.rds"))
coldata <- readRDS(file.path(DIR_PROCESSED, "coldata.rds"))
vsd <- vsd[, rownames(coldata), drop = FALSE]          # alinha colunas->coldata
y   <- factor(coldata$genotype, levels = c("WT", "KO"))

TOPK <- 50   # nº de DEGs selecionados por fold (mesmo tamanho do baseline original)

# --- Seleção de DEGs por fold (limma), controlando o estímulo -------------
# Roda só nas amostras de treino -> sem vazamento da amostra de teste.
select_degs <- function(train_idx, k = TOPK) {
  cd     <- droplevels(coldata[train_idx, ])
  # Se o fold deixar o estímulo sem contraste, cai p/ design só com genótipo.
  design <- tryCatch(model.matrix(~ stimulation + genotype, data = cd),
                     error = function(e) model.matrix(~ genotype, data = cd))
  if (!"genotypeKO" %in% colnames(design))
    design <- model.matrix(~ genotype, data = cd)
  fit <- eBayes(lmFit(vsd[, train_idx, drop = FALSE], design))
  tt  <- topTable(fit, coef = "genotypeKO", number = k, sort.by = "p")
  rownames(tt)
}

finalize_metrics <- function(probs, y) {
  pred <- factor(ifelse(probs >= 0.5, "KO", "WT"), levels = c("WT", "KO"))
  cm   <- caret::confusionMatrix(pred, y, positive = "KO")
  roc  <- pROC::roc(response = y, predictor = probs, levels = c("WT", "KO"), quiet = TRUE)
  list(auc = as.numeric(pROC::auc(roc)),
       precision = cm$byClass["Precision"], recall = cm$byClass["Recall"],
       f1 = cm$byClass["F1"], confusion = cm$table, roc = roc)
}

# --- (A) LOO com seleção de DEGs ANINHADA (estimativa honesta) -------------
loo_nested <- function(y) {
  n <- ncol(vsd); probs <- numeric(n)
  for (i in seq_len(n)) {
    genes_i <- select_degs(setdiff(seq_len(n), i))
    Xtr <- t(vsd[genes_i, -i, drop = FALSE])
    Xte <- t(vsd[genes_i,  i, drop = FALSE])
    fit <- cv.glmnet(Xtr, y[-i], family = "binomial",
                     alpha = 0.5, nfolds = min(5, n - 1))
    probs[i] <- predict(fit, Xte, s = "lambda.min", type = "response")
  }
  finalize_metrics(probs, y)
}

# --- (B) LOO com conjunto FIXO de atributos (eigengenes, não supervisionado)
loo_fixed <- function(X, y) {
  n <- nrow(X); probs <- numeric(n)
  for (i in seq_len(n)) {
    fit <- cv.glmnet(X[-i, , drop = FALSE], y[-i], family = "binomial",
                     alpha = 0.5, nfolds = min(5, n - 1))
    probs[i] <- predict(fit, X[i, , drop = FALSE], s = "lambda.min", type = "response")
  }
  finalize_metrics(probs, y)
}

# Conjunto B: eigengenes dos módulos (rede)
wg <- readRDS(file.path(DIR_PROCESSED, "wgcna.rds"))
Xb <- as.matrix(wg$MEs[rownames(coldata), , drop = FALSE])

res_A <- suppressWarnings(loo_nested(y))      # baseline DEGs — seleção aninhada
res_B <- suppressWarnings(loo_fixed(Xb, y))   # eigengenes WGCNA

metrics <- data.frame(
  modelo   = c("Baseline (DEGs, seleção aninhada)", "Rede (eigengenes WGCNA)"),
  AUC      = c(res_A$auc, res_B$auc),
  Precisao = c(res_A$precision, res_B$precision),
  Recall   = c(res_A$recall, res_B$recall),
  F1       = c(res_A$f1, res_B$f1)
)
write.csv(metrics, file.path(DIR_TABLES, "classificacao_metricas.csv"), row.names = FALSE)
print(metrics)

capture.output(
  cat("== Baseline (DEGs, seleção aninhada no LOO) ==\n"), print(res_A$confusion),
  cat("\n== Rede (eigengenes) ==\n"), print(res_B$confusion),
  file = file.path(DIR_TABLES, "matriz_confusao.txt")
)

png(file.path(DIR_FIGURES, "roc_comparativo.png"), width = 1100, height = 1000, res = 150)
plot(res_A$roc, col = "grey40", main = "ROC — APTX-KO x WT (LOO-CV)")
plot(res_B$roc, col = "firebrick", add = TRUE)
legend("bottomright",
       legend = c(sprintf("Baseline DEGs (AUC=%.2f)", res_A$auc),
                  sprintf("Rede eigengenes (AUC=%.2f)", res_B$auc)),
       col = c("grey40", "firebrick"), lwd = 2)
dev.off()

message("Classificação concluída (seleção de DEGs aninhada no LOO). ",
        "Métricas em results/tables/classificacao_metricas.csv")
