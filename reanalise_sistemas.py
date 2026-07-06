#!/usr/bin/env python3
"""
Reanalise APTX/AOA1 -- metodos de Biologia de Sistemas
Dados publicos: GEO GSE245766 (HMC3 microglia, APTX-KO vs WT x estimulo imune)

Reproduz:
  (1) Discretizacao (quantizacao binaria e ternaria, z-score)
  (2) Estado booleano de consenso por genotipo + distancia de Hamming
  (3) Inferencia de rede regulatoria: correlacao de Pearson vs informacao mutua
  (4) Validacao contra rede de referencia STRING (AUC, AUPR, precisao, recall, F1)

Requer: numpy pandas scipy scikit-learn matplotlib ; acesso a string-db.org
Uso:    python reanalise_sistemas.py GSE245766_counts_table_raw.csv.gz
"""
import sys, urllib.request, urllib.parse
import numpy as np, pandas as pd
from itertools import combinations
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             precision_recall_fscore_support)
from sklearn.feature_selection import mutual_info_regression

RAW = sys.argv[1] if len(sys.argv) > 1 else "GSE245766_counts_table_raw.csv.gz"

# ---- painel de 30 genes: reparo de DNA + imunidade inata/interferon ----
PANEL = {
 "ENSG00000137074":"APTX","ENSG00000042088":"TDP1","ENSG00000073050":"XRCC1",
 "ENSG00000152422":"XRCC4","ENSG00000005156":"LIG3","ENSG00000174405":"LIG4",
 "ENSG00000039650":"PNKP","ENSG00000107290":"SETX","ENSG00000143799":"PARP1",
 "ENSG00000149311":"ATM",
 "ENSG00000137965":"IFI44L","ENSG00000137959":"IFI44","ENSG00000126709":"IFI6",
 "ENSG00000185745":"IFIT1","ENSG00000119917":"IFIT3","ENSG00000089127":"OAS1",
 "ENSG00000111335":"OAS2","ENSG00000135114":"OASL","ENSG00000187608":"ISG15",
 "ENSG00000157601":"MX1","ENSG00000183486":"MX2","ENSG00000134321":"RSAD2",
 "ENSG00000135899":"SP110","ENSG00000185201":"IFITM1","ENSG00000130303":"BST2",
 "ENSG00000115267":"IFIH1","ENSG00000107201":"RIGI","ENSG00000115415":"STAT1",
 "ENSG00000170581":"STAT2","ENSG00000185507":"IRF7",
}

# ---- normalizacao log-CPM ----
raw = pd.read_csv(RAW, index_col=0)
samples = list(raw.columns)
cpm = raw.div(raw.sum(0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)
present = [g for g in PANEL if g in raw.index]
X = logcpm.loc[present].copy(); X.index = [PANEL[g] for g in present]; X = X[samples]

# ---- (1) discretizacao ----
Z = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
B = (Z > 0).astype(int)                                   # binaria
T = pd.DataFrame(np.select([Z > .5, Z < -.5], [1, -1], 0),  # ternaria
                 index=Z.index, columns=Z.columns)

# ---- (2) estado booleano de consenso + Hamming ----
ko = [c for c in samples if c.startswith("APTX_KO")]
wt = [c for c in samples if c.startswith("Control")]
cons_ko = (B[ko].mean(1) > .5).astype(int)
cons_wt = (B[wt].mean(1) > .5).astype(int)
flipped = list(cons_ko.index[cons_ko != cons_wt])
print(f"[Hamming] KO vs WT = {len(flipped)}/{len(present)}  genes: {flipped}")

# ---- STRING gold standard ----
syms = list(X.index)
base = "https://string-db.org/api"
q = urllib.parse.urlencode({"identifiers": "%0d".join(syms), "species": 9606})
net = urllib.request.urlopen(f"{base}/tsv/network?{q}", timeout=60).read().decode()
rows = [l.split("\t") for l in net.strip().split("\n")]
h = rows[0]; ai, bi, si = h.index("preferredName_A"), h.index("preferredName_B"), h.index("score")
gold = {frozenset((r[ai], r[bi])) for r in rows[1:]
        if r[ai] in syms and r[bi] in syms and float(r[si]) >= 0.4}
print(f"[STRING] arestas gold (score>=0.4) = {len(gold)}")

# ---- (3) inferencia: correlacao vs informacao mutua ----
Xv = X.values; n = len(syms)
pairs = list(combinations(range(n), 2))
y = np.array([1 if frozenset((syms[i], syms[j])) in gold else 0 for i, j in pairs])

C = np.corrcoef(Xv)
score_corr = np.array([abs(C[i, j]) for i, j in pairs])

MI = np.zeros((n, n))
for t in range(n):
    MI[t] = mutual_info_regression(Xv.T, Xv[t], random_state=0)
MI = (MI + MI.T) / 2
score_mi = np.array([MI[i, j] for i, j in pairs])

# ---- (4) metricas ----
def pt(scores, k=None):
    k = int(y.sum()) if k is None else k
    pred = np.zeros_like(y); pred[np.argsort(-scores)[:k]] = 1
    p, r, f1, _ = precision_recall_fscore_support(y, pred, average="binary", zero_division=0)
    return p, r, f1

for name, sc in [("Correlacao(|r|)", score_corr), ("Info.Mutua", score_mi)]:
    p, r, f1 = pt(sc)
    print(f"[{name:16s}] AUC={roc_auc_score(y,sc):.3f} "
          f"AUPR={average_precision_score(y,sc):.3f} P={p:.3f} R={r:.3f} F1={f1:.3f}")
print(f"[Aleatorio]        AUPR(prevalencia)={y.mean():.3f}  pares={len(y)}")
