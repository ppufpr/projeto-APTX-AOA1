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
R/                          Análises em R/Bioconductor (DESeq2, WGCNA, CoDiNA) — ver R/README.md
grn_metrics.json            Métricas de validação (correlação vs informação mútua)
figure_reanalysis.png       Fig. 1: DEGs, escore ISG, teste de permutação
figure_discretizacao.png    Fig. 2: quantização ternária + estado booleano (Hamming)
figure_grn_inference.png    Fig. 3: ROC/PR + métricas (inferência vs STRING)
```
O diretório `R/` reúne o código em R que fez parte do desenvolvimento do artigo:
as etapas de expressão diferencial (**DESeq2**), co-expressão (**WGCNA**) e
co-expressão diferencial (**CoDiNA**) descritas nos Resultados. As demais etapas
(discretização, estado booleano, inferência de rede e validação STRING) estão em
`reanalise_sistemas.py`.

## Reprodução
```bash
# dados (uma vez):
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245766/suppl/GSE245766_counts_table_raw.csv.gz

pip install numpy pandas scipy scikit-learn matplotlib
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
pdflatex artigo_APTX_AOA1.tex   # 2x para referências cruzadas
# ou: tectonic artigo_APTX_AOA1.tex
```

## Informação mútua
Informação mútua é uma medida estatística que indica quanto conhecer uma variável ajuda a reduzir a incerteza sobre outra variável.

## Fonte dos dados
Madsen HB, et al. *The DNA repair enzyme, Aprataxin, plays a role in innate
immune signaling.* Front Aging Neurosci. 2023;15:1290681.
doi:10.3389/fnagi.2023.1290681. Dados: GEO GSE245766.
