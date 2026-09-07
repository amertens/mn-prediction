"""
scripts/protocol_v2/51_hapc_smoke_test.py   [HP-01]

Smoke test of the hapc package (Wang, Schuler, van der Laan, Garcia Meixide,
arXiv 2602.10613: principal-component Highly Adaptive Lasso / Ridge) on the
protocol's own design matrices, scored the way the protocol scores its arms.

  (a) in-fill, Ghana child iron, 75 districts x 97 domain PCs: honest 5-fold
      out-of-fold predictions under the protocol's own district folds (three
      of the ten draws), lambda tuned by hapc's inner CV inside each training
      fold; norms "2" (PCHAR, closed form), "1" (PCHAL, soft-threshold) and
      "sv" (projected gradient descent on the sectional-variation ball).
  (b) leave-one-country-out, child iron and child vitamin A, 206 districts x
      123 domain PCs: fit on three countries, predict the fourth.
Reference: the protocol's zero-tuning domain index on the same cells
(benchmarks_v2_cells.csv), and a plain ridge on the same PCs as a floor.

Run with the reticulate venv python (hapc 2.6.0 installed there):
  ".virtualenvs/r-reticulate/Scripts/python.exe" scripts/protocol_v2/51_hapc_smoke_test.py
-> results/tables/protocol_v2/hapc_smoke/hapc_smoke_results.csv
"""
import os, sys, time, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
warnings.filterwarnings("ignore")
import hapc
from hapc.cv import cv_hapc

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
D = ROOT + "results/tables/protocol_v2/hapc_smoke/"
MAXDEG = int(os.environ.get("HAPC_MAXDEG", "1"))
NORMS = os.environ.get("HAPC_NORMS", "2,1,sv").split(",")
REPS = int(os.environ.get("HAPC_REPS", "3"))
GRID = dict(log_lambda_min=float(os.environ.get("HAPC_LMIN", "-8")), log_lambda_max=float(os.environ.get("HAPC_LMAX", "2")), grid_length=int(os.environ.get("HAPC_GRIDN", "12")))
rows = []

def pc_cols(df):
    return [c for c in df.columns if c not in ("country", "Admin1", "Admin2", "y_mod", "y_nat", "w") and not c.startswith("fold_rep")]

def fit_predict(Xtr, ytr, Xte, norm):
    """hapc with inner CV for lambda on the training rows; returns predictions for Xte and the chosen lambda."""
    t0 = time.time()
    res = cv_hapc(Xtr, ytr, family="gaussian", max_degree=MAXDEG, npcs=None, norm=norm, nfolds=5, predict=Xte, **GRID)
    pred = np.asarray(res.predictions).ravel()
    lam = getattr(res, "best_lambda", np.nan)
    return pred, lam, time.time() - t0

def ridge_predict(Xtr, ytr, Xte):
    m = RidgeCV(alphas=np.logspace(-3, 4, 30)).fit(Xtr, ytr)
    return m.predict(Xte)

def rho(y, p):
    if np.nanstd(p) == 0: return 0.0
    return spearmanr(y, p, nan_policy="omit").correlation

# ── (a) in-fill ─────────────────────────────────────────────────────────────
for target in ("level", "prev"):
    df = pd.read_csv(D + f"infill_ghana_child_iron_{target}.csv")
    X = df[pc_cols(df)].to_numpy(float); y = df["y_mod"].to_numpy(float); ynat = df["y_nat"].to_numpy(float)
    for norm in NORMS + ["ridge_pcs"]:
        rhos, lams, secs = [], [], []
        for rep in range(1, REPS + 1):
            folds = df[f"fold_rep{rep}"].to_numpy(int); pred = np.full(len(y), np.nan)
            for f in np.unique(folds):
                te = folds == f; tr = ~te
                if norm == "ridge_pcs":
                    pred[te] = ridge_predict(X[tr], y[tr], X[te]); lam, sec = np.nan, 0.0
                else:
                    pred[te], lam, sec = fit_predict(X[tr], y[tr], X[te], norm)
                lams.append(lam); secs.append(sec)
            rhos.append(rho(ynat, pred))
        rows.append(dict(estimand="infill", cell=f"Ghana child_iron", target=target, model=norm, max_degree=MAXDEG, n=len(y), p=X.shape[1],
                         spearman_median=float(np.median(rhos)), spearman_min=float(np.min(rhos)), spearman_max=float(np.max(rhos)),
                         median_lambda=float(np.nanmedian(lams)) if len(lams) else np.nan, seconds_per_fit=float(np.mean(secs))))
        print(f"infill {target:5s} {norm:9s} rho median {np.median(rhos):+.3f} (draws {np.min(rhos):+.3f}..{np.max(rhos):+.3f}) | {np.mean(secs):.1f} s per fit", flush=True)

# ── (b) leave-one-country-out ───────────────────────────────────────────────
for on in ("child_iron", "child_vitA"):
    for target in ("level", "prev"):
        df = pd.read_csv(D + f"loco_{on}_{target}.csv")
        X = df[pc_cols(df)].to_numpy(float); y = df["y_mod"].to_numpy(float); ynat = df["y_nat"].to_numpy(float); ctry = df["country"].to_numpy()
        for norm in NORMS + ["ridge_pcs"]:
            out = {}
            for cn in np.unique(ctry):
                te = ctry == cn; tr = ~te
                if norm == "ridge_pcs":
                    pred = ridge_predict(X[tr], y[tr], X[te]); lam, sec = np.nan, 0.0
                else:
                    pred, lam, sec = fit_predict(X[tr], y[tr], X[te], norm)
                out[cn] = rho(ynat[te], pred)
                rows.append(dict(estimand="country", cell=f"{cn} {on}", target=target, model=norm, max_degree=MAXDEG, n=int(te.sum()), p=X.shape[1],
                                 spearman_median=float(out[cn]), spearman_min=np.nan, spearman_max=np.nan, median_lambda=float(lam), seconds_per_fit=float(sec)))
            print(f"LOCO {on:10s} {target:5s} {norm:9s} " + " | ".join(f"{k} {v:+.2f}" for k, v in out.items()) + f" | mean {np.mean(list(out.values())):+.3f}", flush=True)

pd.DataFrame(rows).to_csv(D + "hapc_smoke_results.csv", index=False)
print("DONE", flush=True)
