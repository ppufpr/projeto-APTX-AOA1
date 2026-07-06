# ============================================================
# 08_network_medicine.R — Módulo das ataxias de reparo de DNA (STRING)
# Constrói a rede de interações (STRINGdb) dos genes-semente
# (APTX, SETX, PNKP, ATM, TDP1, XRCC1, XRCC4, ...), mede centralidade
# e verifica se o APTX é um hub conectado ao módulo de reparo de DNA.
# Integra com a expressão: marca quais genes do módulo são DEGs no KO.
# ============================================================
source("R/config.R")
suppressMessages({ library(STRINGdb); library(igraph) })

# 1. Inicializa STRINGdb (humano)
string_db <- STRINGdb$new(version = "12.0", species = STRING_TAX,
                          score_threshold = STRING_SCORE,
                          input_directory = DIR_RAW)

# 2. Mapeia genes-semente para IDs do STRING
seed_df <- data.frame(gene = SEED_GENES)
mapped  <- string_db$map(seed_df, "gene", removeUnmappedRows = TRUE)
message("Genes-semente mapeados no STRING: ", nrow(mapped), "/", length(SEED_GENES))

# 3. Vizinhança: interatores diretos dos genes-semente
neighbors <- string_db$get_neighbors(mapped$STRING_id)
all_ids   <- unique(c(mapped$STRING_id, neighbors))

# 4. Sub-rede + igraph
edges <- string_db$get_interactions(all_ids)
if (nrow(edges) > 0) {
  g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE)
  g <- simplify(g)

  # Centralidades
  cent <- data.frame(
    node       = V(g)$name,
    degree     = degree(g),
    betweenness= betweenness(g),
    closeness  = closeness(g)
  )
  # Marca os genes-semente
  cent$is_seed <- cent$node %in% mapped$STRING_id
  cent <- cent[order(-cent$degree), ]
  write.csv(cent, file.path(DIR_TABLES, "netmed_centrality.csv"), row.names = FALSE)

  # 5. Integração com expressão: quais nós são alterados no KO?
  #    Dois critérios: (a) DEG estrito (|log2FC|>1 & padj<0.05) e
  #    (b) significativo por padj (qualquer LFC) — os genes de reparo mudam
  #    pouco em magnitude, mas vários são significativos (APTX, TDP1, XRCC1...).
  deg_file <- file.path(DIR_TABLES, "de_genotype_KO_vs_WT_degs.csv")
  all_file <- file.path(DIR_TABLES, "de_genotype_KO_vs_WT_all.csv")
  if (file.exists(deg_file)) {
    degs <- read.csv(deg_file)$gene
    deg_map <- string_db$map(data.frame(gene = degs), "gene", removeUnmappedRows = TRUE)
    cent$is_DEG_KO <- cent$node %in% deg_map$STRING_id

    if (file.exists(all_file)) {
      allres <- read.csv(all_file)
      idcol  <- if ("gene" %in% names(allres)) "gene" else names(allres)[1]
      signif <- allres[[idcol]][!is.na(allres$padj) & allres$padj < DE_PADJ]
      sig_map <- string_db$map(data.frame(gene = signif), "gene", removeUnmappedRows = TRUE)
      cent$is_signif_KO <- cent$node %in% sig_map$STRING_id
    }
    write.csv(cent, file.path(DIR_TABLES, "netmed_centrality.csv"), row.names = FALSE)
    message("Genes-semente que são DEGs estritos (|LFC|>1 & padj<", DE_PADJ, "): ",
            sum(cent$is_seed & cent$is_DEG_KO))
    if ("is_signif_KO" %in% names(cent))
      message("Genes-semente significativos por padj<", DE_PADJ, " (qualquer LFC): ",
              sum(cent$is_seed & cent$is_signif_KO))
  }

  # 6. Figura da sub-rede dos genes-semente
  png(file.path(DIR_FIGURES, "netmed_seed_module.png"), width = 1200, height = 1200, res = 140)
  string_db$plot_network(mapped$STRING_id)
  dev.off()

  message("Medicina de Redes concluída. Centralidades em results/tables/netmed_centrality.csv")
} else {
  message("Nenhuma aresta retornada — verifique conexão ao STRING e o STRING_SCORE.")
}
