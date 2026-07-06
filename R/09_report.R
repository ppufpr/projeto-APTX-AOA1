# ============================================================
# 09_report.R — Gera o relatório (.docx) com figuras, tabelas e resultados.
# Usa officer (Word nativo, sem pandoc) + flextable.
# Saída: results/relatorio_APTX_AOA1.docx
# IMPORTANTE: officer rastreia o cursor no objeto retornado -> SEMPRE reatribuir
# (aqui via helpers com `doc <<- ...`), senão o documento sai em ordem invertida.
# ============================================================
source("R/config.R")
suppressMessages({
  library(officer); library(flextable)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

rd  <- function(f) read.csv(file.path(DIR_TABLES, f))
sym <- function(ens) { s <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL",
                                    column = "SYMBOL", multiVals = "first"); ifelse(is.na(s), ens, s) }
pe  <- function(x) formatC(x, format = "e", digits = 2)         # p-valor científico
FIG <- function(f) file.path(DIR_FIGURES, f)

png_ratio <- function(f) { con <- file(f, "rb"); on.exit(close(con)); readBin(con, "raw", 16)
  w <- readBin(con, "integer", 1, 4, endian = "big"); h <- readBin(con, "integer", 1, 4, endian = "big"); h / w }
ft_std <- function(df) { f <- flextable(df); f <- theme_booktabs(f); f <- fontsize(f, size = 9, part = "all")
  f <- autofit(f); f <- set_table_properties(f, layout = "autofit", width = 1); f }

doc <- read_docx()

# Helpers que REATRIBUEM doc no escopo global (avançam o cursor corretamente)
H1  <- function(t) doc <<- body_add_par(doc, t, style = "heading 1")
H2  <- function(t) doc <<- body_add_par(doc, t, style = "heading 2")
P   <- function(t) doc <<- body_add_par(doc, t, style = "Normal")
CAP <- function(t) doc <<- body_add_par(doc, t, style = "centered")
FT  <- function(df) doc <<- body_add_flextable(doc, ft_std(df))
BRK <- function() doc <<- body_add_break(doc)
FIGURE <- function(file, w = 6) doc <<- body_add_img(doc, FIG(file), width = w,
                                                     height = round(w * png_ratio(FIG(file)), 2))
FPAR <- function(txt, size, bold = FALSE, italic = FALSE)
  doc <<- body_add_fpar(doc, fpar(ftext(txt, fp_text(bold = bold, italic = italic, font.size = size)),
                                  fp_p = fp_par(text.align = "center")))

# ============================================================
# CAPA
# ============================================================
FPAR("Redes gênicas e a perda de APTX (AOA1):", 20, bold = TRUE)
FPAR("reprogramação do transcriptoma e o módulo das ataxias de reparo de DNA", 14, italic = TRUE)
P("")
CAP("Projeto final — Biologia de Sistemas (PPGAB / UFPR–UTFPR)")
CAP("Dados: GSE245766 — microglia humana HMC3, desenho fatorial 2×2 (genótipo APTX-KO/WT × estímulo imune IS/NS), n = 12")
CAP(format(Sys.Date(), "%d/%m/%Y"))
BRK()

# ============================================================
# RESUMO
# ============================================================
H1("Resumo")
P(paste0(
  "A aprataxina (APTX) repara quebras de fita simples de DNA e é a causa da Ataxia com Apraxia ",
  "Oculomotora tipo 1 (AOA1). Investigamos como a perda de APTX reorganiza o transcriptoma de microglia ",
  "humana (HMC3) num desenho fatorial 2×2 (genótipo × estímulo imune; n = 12), por expressão ",
  "diferencial (DESeq2), co-expressão (WGCNA), co-expressão diferencial (CoDiNA), enriquecimento ",
  "funcional (GO/KEGG), um classificador validado (LOO-CV com seleção de atributos aninhada) e Medicina ",
  "de Redes (STRINGdb). O efeito do genótipo abrange 1.221 genes diferencialmente expressos (DEGs), ",
  "dominados pelo programa de interferon/antiviral (IFI44L, IFITM1, OAS1, RIG-I/DDX58), majoritariamente ",
  "reprimido no KO. O enriquecimento destaca resposta a vírus, sinalização por interferon tipo I/II, ",
  "NF-κB e receptores NOD-like. Os próprios genes de reparo de DNA mudam pouco em magnitude (nenhum ",
  "cruza |log2FC| > 1), mas vários são estatisticamente significativos. Concluímos que a perda de APTX ",
  "não reprograma transcricionalmente os genes de reparo em si — atua no nível enzimático —, mas ",
  "reconfigura a jusante o programa imune inato, coerente com o papel emergente do APTX na via ",
  "cGAS-STING / RIG-I-MAVS."))

# ============================================================
# 1. INTRODUÇÃO
# ============================================================
H1("1. Introdução")
P(paste0(
  "O APTX (aprataxina) resolve intermediários abortivos 5′-adenilados da ligação de DNA, atuando no ",
  "reparo de quebras de fita simples (SSBR) em parceria com XRCC1/XRCC4. Sua perda causa a AOA1, ataxia ",
  "recessiva do grupo das ataxias de reparo de DNA (com SETX/AOA2, PNKP/AOA4, ATM/AT). Um papel ",
  "recém-descrito na imunidade inata (cGAS-STING, RIG-I/MAVS) torna o gene um alvo rico para análise ",
  "de redes. Pergunta: como a perda de APTX reorganiza a rede de co-expressão gênica e quais processos ",
  "(reparo de DNA, resposta imune inata) são afetados? O APTX e seus parceiros formam um módulo coerente ",
  "na rede de interações?"))

# ============================================================
# 2. DADOS E MÉTODOS
# ============================================================
H1("2. Dados e Métodos")
P(paste0(
  "Dados: GSE245766 (RNA-seq), linhagem de microglia humana HMC3, desenho fatorial 2×2 — genótipo ",
  "(APTX-KO vs. selvagem) × estímulo imune (estimulado IS vs. não-estimulado NS), 3 réplicas por ",
  "condição (n = 12). Pré-processamento: filtro de baixa expressão (edgeR::filterByExpr; 15.600 genes ",
  "retidos) e estabilização de variância (DESeq2::vst)."))
P(paste0(
  "Expressão diferencial: DESeq2 com desenho ~ genótipo + estímulo + genótipo:estímulo (corte: ",
  "padj < 0,05 e |log2FC| > 1). Co-expressão: WGCNA (rede signed). Como o ajuste scale-free falha com ",
  "n pequeno, usou-se o power recomendado pelo WGCNA para redes signed (β = 18). Co-expressão ",
  "diferencial: CoDiNA (rede APTX-KO × WT). Enriquecimento: clusterProfiler (GO BP e KEGG, humano). ",
  "Classificação KO×WT: glmnet (elastic-net) com validação leave-one-out; para evitar viés de ",
  "seleção (Ambroise & McLachlan, 2002), a seleção de DEGs é refeita dentro de cada fold (limma), ",
  "enquanto os eigengenes do WGCNA — não supervisionados — entram como conjunto fixo. Medicina de ",
  "Redes: STRINGdb v12.0 (Homo sapiens, score ≥ 400) + igraph (grau, intermediação, proximidade)."))

# ============================================================
# 3. RESULTADOS
# ============================================================
H1("3. Resultados")

## 3.1 Expressão diferencial
H2("3.1 Expressão diferencial (DESeq2)")
deg_counts <- data.frame(
  Contraste = c("Genótipo (KO vs WT)", "Estímulo (IS vs NS)", "Interação genótipo×estímulo"),
  Testados  = c(15600, 15600, 15600),
  DEGs      = c(1221, 1291, 209),
  Up        = c(554, 936, 82),
  Down      = c(667, 355, 127))
FT(deg_counts)
CAP("Tabela 1. Número de genes diferencialmente expressos por contraste (padj < 0,05 e |log2FC| > 1).")
P("")
g <- rd("de_genotype_KO_vs_WT_degs.csv"); g <- g[order(g$padj, -abs(g$log2FoldChange)), ]
top12 <- data.frame(Gene = sym(head(g$gene, 12)),
                    log2FC = round(head(g$log2FoldChange, 12), 2),
                    padj = pe(head(g$padj, 12)))
FT(top12)
CAP("Tabela 2. Top 12 DEGs do efeito do genótipo (KO vs WT). Predomina o programa interferon/antiviral.")
P("")
FIGURE("volcano_genotype.png", w = 5.5)
CAP("Figura 1. Volcano plot do efeito do genótipo (APTX-KO vs WT).")
P(paste0(
  "O efeito do genótipo abrange 1.221 DEGs. Entre os mais significativos estão genes estimulados por ",
  "interferon (IFI44L, IFITM1, OAS1, IFI6, IFIT1, USP18, BST2) e o sensor citosólico RIG-I (DDX58), ",
  "majoritariamente reprimidos no KO — já apontando o eixo imune inato."))

## 3.2 WGCNA
H2("3.2 Co-expressão (WGCNA)")
FIGURE("wgcna_scalefree.png", w = 6.5)
CAP("Figura 2. Ajuste scale-free e conectividade média vs. power. Com n = 12 o ajuste é fraco; adotou-se β = 18 (recomendação WGCNA para redes signed).")
P("")
mt <- rd("wgcna_module_trait.csv")
mt2 <- data.frame(Modulo = sub("^ME", "", mt$module),
                  cor_genotipo = round(mt$cor_genotype, 2), p_genotipo = pe(mt$p_genotype),
                  cor_estimulo = round(mt$cor_stim, 2),    p_estimulo = pe(mt$p_stim))
FT(mt2)
CAP("Tabela 3. Correlação módulo–traço (eigengene vs. genótipo e estímulo).")
P(paste0(
  "Foram detectados 8 módulos. Os módulos yellow (r = -0,98; p = 1,3×10⁻⁸) e blue (r = +0,96; ",
  "p = 9,3×10⁻⁷) associam-se fortemente ao genótipo, enquanto turquoise (r = +0,98) e red (r = -0,96) ",
  "captam o estímulo — ou seja, o desenho fatorial separa-se bem em módulos de genótipo e de estímulo. ",
  "Dos genes-semente, apenas SETX (turquoise), TDP1 e LIG3 (blue) entram na rede; os demais (incl. APTX) ",
  "ficam fora do conjunto dos 5.000 genes mais variáveis — consistente com sua baixa variação transcricional."))

## 3.3 CoDiNA
H2("3.3 Co-expressão diferencial (CoDiNA)")
codina <- data.frame(
  Categoria = c("Comum (α)", "Específico KO (g.KO)", "Específico WT (g.WT)", "Específico WT (b.WT)", "Indefinido (U)"),
  Nos = c(1454, 224, 291, 27, 1))
FT(codina)
CAP("Tabela 4. Classificação dos nós na rede de co-expressão diferencial (CoDiNA), KO vs WT.")
P(paste0(
  "A comparação das redes de co-expressão revela rewiring substancial: 224 nós têm conectividade ",
  "específica do KO e 291 específica do WT, contra 1.454 conservados. Ou seja, a perda de APTX não ",
  "muda apenas níveis de expressão, mas a própria arquitetura de co-expressão."))

## 3.4 Enriquecimento
H2("3.4 Enriquecimento funcional (GO/KEGG)")
go <- rd("enrich_GO_BP.csv")
go12 <- data.frame(Termo_GO_BP = substr(go$Description[1:12], 1, 55),
                   p_ajust = pe(go$p.adjust[1:12]), Genes = go$Count[1:12])
FT(go12)
CAP("Tabela 5. Top 12 termos GO (Processo Biológico) enriquecidos (de 368 termos).")
P("")
FIGURE("enrich_GO_BP.png", w = 6)
CAP("Figura 3. Dotplot dos termos GO BP enriquecidos (APTX-KO vs WT).")
P("")
k <- rd("enrich_KEGG.csv")
k12 <- data.frame(Via_KEGG = substr(k$Description[1:min(12,nrow(k))], 1, 55),
                  p_ajust = pe(k$p.adjust[1:min(12,nrow(k))]), Genes = k$Count[1:min(12,nrow(k))])
FT(k12)
CAP("Tabela 6. Vias KEGG enriquecidas (de 23 vias).")
P("")
FIGURE("enrich_KEGG.png", w = 6)
CAP("Figura 4. Dotplot das vias KEGG enriquecidas.")
P(paste0(
  "GO e KEGG convergem no eixo imune inato/antiviral: resposta a vírus, sinalização por interferon ",
  "tipo I e II, produção de IFN-β, e vias KEGG de NF-κB, receptores NOD-like, Influenza A e Measles. ",
  "É a assinatura cGAS-STING / RIG-I-MAVS esperada."))

## 3.5 Classificação
H2("3.5 Classificação (LOO-CV) e métricas")
cm <- rd("classificacao_metricas.csv")
FT(cm)
CAP("Tabela 7. Métricas de classificação KO vs WT (LOO-CV; seleção de DEGs aninhada).")
P("")
FIGURE("roc_comparativo.png", w = 4.8)
CAP("Figura 5. Curvas ROC — baseline (DEGs) vs. rede (eigengenes WGCNA).")
P(paste0(
  "Ambos os modelos atingem AUC = 1 e classificação perfeita (12/12). Como a seleção de DEGs foi ",
  "aninhada no LOO (sem vazamento), a separação perfeita reflete a forte distinção biológica entre ",
  "WT e knockout, não otimismo metodológico. A contrapartida é que a tarefa é fácil demais para ",
  "discriminar os dois conjuntos de atributos."))

## 3.6 Medicina de Redes
H2("3.6 Medicina de Redes (STRINGdb)")
FIGURE("netmed_seed_module.png", w = 5)
CAP("Figura 6. Sub-rede STRING dos genes-semente (módulo das ataxias de reparo de DNA).")
P("")
cent <- rd("netmed_centrality.csv")
seeds <- c("APTX","SETX","PNKP","ATM","TDP1","XRCC1","XRCC4","LIG4","LIG3","PARP1")
suppressMessages(library(STRINGdb))
sdb <- STRINGdb$new(version = "12.0", species = STRING_TAX, score_threshold = STRING_SCORE, input_directory = DIR_RAW)
smap <- sdb$map(data.frame(gene = seeds), "gene", removeUnmappedRows = TRUE)
sc <- merge(cent[cent$is_seed, ], smap, by.x = "node", by.y = "STRING_id")
sc <- sc[order(-sc$degree), ]
cent_tab <- data.frame(Gene = sc$gene, Grau = sc$degree,
                       Intermediacao = round(sc$betweenness, 0),
                       Signif_KO_padj = ifelse(sc$is_signif_KO, "sim", "não"))
FT(cent_tab)
CAP("Tabela 8. Centralidade dos genes-semente na rede STRING e significância no KO (padj < 0,05).")
P(paste0(
  "ATM e PARP1 são os hubs do módulo (maior grau e intermediação), reguladores centrais do reparo. ",
  "O APTX, embora de grau modesto, exibe intermediação alta (~12.430), atuando como conector/bottleneck ",
  "entre subgrupos — coerente com seu papel de ponte no SSBR. Nenhum gene-semente cruza o corte estrito ",
  "de |log2FC| > 1, mas cinco são significativos por padj (APTX, TDP1, XRCC1, XRCC4, LIG3): o eixo ",
  "funcional ao redor do APTX é sutilmente reprogramado."))

# ============================================================
# 4. DISCUSSÃO
# ============================================================
H1("4. Discussão")
P(paste0(
  "Os resultados convergem para uma narrativa coerente. A perda de APTX não altera substancialmente a ",
  "expressão dos próprios genes de reparo de DNA — esperado, pois o APTX atua no nível enzimático ",
  "(resolução de intermediários de ligação), não como regulador transcricional. O reflexo ",
  "transcricional forte ocorre a jusante, no programa imune inato/interferon, captado de forma ",
  "consistente pela expressão diferencial (RIG-I, OAS1, ISGs), pelo enriquecimento (resposta a vírus, ",
  "IFN tipo I/II, NF-κB, NOD-like) e pelo rewiring da co-expressão (CoDiNA). Isso conecta-se ao papel ",
  "emergente do APTX na via cGAS-STING / RIG-I-MAVS. Na rede de interações, APTX comporta-se como ",
  "conector (alta intermediação) dentro do módulo das ataxias de reparo, ancorado por hubs como ATM e ",
  "PARP1."))

# ============================================================
# 5. LIMITAÇÕES
# ============================================================
H1("5. Limitações")
P(paste0(
  "(i) Modelo celular (microglia HMC3), não tecido de paciente; n pequeno (3/condição), o que torna o ",
  "WGCNA exploratório e justifica o foco em módulos, não arestas isoladas. (ii) O ajuste scale-free ",
  "falha com n = 12; o power foi fixado por recomendação (β = 18), não por estimativa de dados. ",
  "(iii) A classificação KO×WT é trivial (AUC = 1) e não discrimina conjuntos de atributos; uma tarefa ",
  "mais difícil (efeito do estímulo ou interação) seria mais informativa. (iv) O APTX aparece ligeiramente ",
  "aumentado em mRNA (log2FC = +0,27), sugerindo knockout funcional/proteico com mRNA preservado — a ",
  "estratégia de KO do estudo original deve ser confirmada."))

# ============================================================
# 6. CONCLUSÃO
# ============================================================
H1("6. Conclusão")
P(paste0(
  "A perda de APTX em microglia humana reprograma o transcriptoma sobretudo no eixo imune inato/interferon, ",
  "e não pela alteração transcricional dos genes de reparo em si. Os parceiros de APTX no SSBR/NHEJ ",
  "(XRCC1, XRCC4, LIG3, TDP1) mostram mudanças sutis porém significativas, e o APTX ocupa posição de ",
  "conector no módulo das ataxias de reparo de DNA. Os achados situam a AOA1 na interface entre reparo de ",
  "DNA e imunidade inata."))

# Reprodutibilidade
H1("Reprodutibilidade")
P(paste0("Seed fixa (", SEED, "). Pipeline em R/00–09. Ambiente: R ", as.character(getRversion()),
         " + Bioconductor 3.18 (versões fixadas). Dados brutos não versionados (termos GEO/STRING)."))

# ============================================================
# SALVAR
# ============================================================
out <- file.path("results", "relatorio_APTX_AOA1.docx")
print(doc, target = out)
message("Relatório gerado: ", normalizePath(out))
