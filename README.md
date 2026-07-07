# Redes gênicas e a perda de APTX (AOA1) — projeto final de Bioinformática

Reanálise dos dados públicos **GEO GSE245766** (microglia humana HMC3;
APTX-KO vs WT x estímulo imune) aplicando métodos de bioinformática:
discretização, modelagem booleana e inferência de rede regulatória validada
contra rede de referência (STRING).

## Conteúdo
```
artigo_APTX_AOA1.tex        Fonte LaTeX do artigo (compila com pdflatex/tectonic)
artigo_APTX_AOA1.pdf        Artigo final
reanalise_sistemas.py       Script reproduzível (discretização, booleano, GRN, métricas)
R/                          Análises em R/Bioconductor (DESeq2, WGCNA, CoDiNA)
grn_metrics.json            Métricas de validação (correlação vs informação mútua)
figure_reanalysis.png       Fig. 1: DEGs, escore ISG, teste de permutação
figure_discretizacao.png    Fig. 2: quantização ternária + estado booleano (Hamming)
figure_grn_inference.png    Fig. 3: ROC/PR + métricas (inferência vs STRING)
```
### Diretório `R/`
Pipeline em R/Bioconductor que gerou as análises principais do artigo (delineamento
fatorial 2×2, n=12). Rodar `source("R/00_setup.R")` uma vez e depois
`source("R/run_all.R")`:
```
config.R                 Parâmetros centrais (GSE245766, org. humano, táxon STRING)
00_setup.R               Instala/carrega dependências (CRAN + Bioconductor)
01_download.R            Baixa dados e metadados do GEO (GEOquery)
02_preprocess.R          Metadados 2×2, filtro e normalização (edgeR/DESeq2, vst)
03_diff_expression.R     Expressão diferencial fatorial (DESeq2): KO, estímulo, interação
04_wgcna.R               Rede de co-expressão e módulos (WGCNA)
05_diff_coexpression.R   Co-expressão diferencial APTX-KO vs WT (CoDiNA)
06_enrichment.R          Enriquecimento funcional GO/KEGG (clusterProfiler)
07_classification.R      Classificador KO×WT com validação LOO (DEGs vs eigengenes)
08_network_medicine.R    Módulo de reparo de DNA via STRINGdb + centralidade (igraph)
09_report.R              Relatório .docx com figuras e tabelas (officer/flextable)
run_all.R                Executa os passos 01–09 na ordem
```
As demais etapas (discretização, estado booleano, inferência de rede correlação vs.
informação mútua e validação STRING) estão no script Python `reanalise_sistemas.py`.

## Reprodução
```bash
# dados (uma vez):
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245766/suppl/GSE245766_counts_table_raw.csv.gz

## Ambiente

Crie um ambiente venv

pip install numpy pandas scipy scikit-learn matplotlib ou, se preferir, pip install -r requirements.txt
python reanalise_sistemas.py GSE245766_counts_table_raw.csv.gz
```
Requer acesso a `string-db.org` para a rede de referência.

## Resultados principais
| Método | AUC | AUPR | Precisão | Recall | F1 |
|---|---|---|---|---|---|
| Correlação (\|r\|) | 0,70 | 0,71 | 0,63 | 0,63 | 0,63 |
| Informação mútua  | 0,59 | 0,55 | 0,59 | 0,59 | 0,59 |
| Aleatório         | 0,50 | 0,49 | — | — | — |

- Estado booleano de consenso: **distância de Hamming = 10/30** entre WT e KO.
- No cenário de baixa disponibilidade amostral (M = 12), a correlação apresenta desempenho superior ao da informação mútua.

## Compilar o artigo
```bash
pdflatex artigo_APTX_AOA1.tex
```

## Sempre aprendendo
Informação mútua é uma medida estatística que indica quanto conhecer uma variável ajuda a reduzir a incerteza sobre outra variável.

Interferons (IFNs) são proteínas de sinalização (citocinas) que as células liberam quando detectam uma ameaça

## Diretório "illustrations"
Algumas imagens criadas durante a compreensão dos conhecimentos vistos nas bibliográfias encontradas durante a revisão da literatura de conceitos biológicos.

## Fonte dos dados
Madsen HB, et al. *The DNA repair enzyme, Aprataxin, plays a role in innate
immune signaling.* Front Aging Neurosci. 2023;15:1290681.
doi:10.3389/fnagi.2023.1290681. Dados: GEO GSE245766.
