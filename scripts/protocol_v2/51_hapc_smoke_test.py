"""
scripts/protocol_v2/51_hapc_smoke_test.py   [HP-01, HP-02]

Smoke test of the hapc package (Wang, Schuler, van der Laan, Garcia Meixide,
arXiv 2602.10613: principal-component Highly Adaptive Lasso / Ridge) on the
protocol's own design matrices, scored the way the protocol scores its arms.

  (a) in-fill, Ghana child iron, 75 districts: honest 5-fold out-of-fold
      predictions under the protocol's own district folds (three of the ten
      draws), lambda tuned by hapc's inner CV inside each training fold;
      norms "2" (PCHAR, closed form) and "1" (PCHAL, soft-threshold); "sv"
      (projected gradient descent on the sectional-variation ball) on request.
  (b) leave-one-country-out, child iron and child vitamin A, ~206 districts:
      fit on three countries, predict the fourth.

HP-02 (2026-09-17): every cell is run on TWO designs from script 50 --
  pcs  the protocol's domain PCs (`<cell>.csv`; the HP-01 design)
  raw  the rank-normalised predictor matrix itself (`<cell>_raw.csv`), so
       hapc's own kernel-PC reduction competes with the domain PCA.
Reference: the protocol's zero-tuning domain index on the same cells
(benchmarks_v2_cells.csv), and a plain ridge on the same design as a floor.

Run with the reticulate venv python (hapc 2.6.0 installed there):
  ".virtualenvs/r-reticulate/Scripts/python.exe" scripts/protocol_v2/51_hapc_smoke_test.py
Env: HAPC_DESIGNS (pcs,raw), HAPC_NORMS (2,1), HAPC_REPS (3), HAPC_MAXDEG (1),
     HAPC_LMIN/LMAX/GRIDN (lambda grid), HAPC_OUT (hapc_design_comparison.csv)
-> results/tables/protocol_v2/hapc_smoke/hapc_design_comparison.csv
   (HP-01's hapc_smoke_results.csv is left as the record of that run)
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
NORMS = os.environ.get("HAPC_NORMS", "2,1").split(",")
DESIGNS = os.environ.get("HAPC_DESIGNS", "pcs,raw").split(",")
REPS = int(os.environ.get("HAPC_REPS", "3"))
OUT = os.environ.get("HAPC_OUT", "hapc_design_comparison.csv")
GRID = dict(log_lambda_min=float(os.environ.get("HAPC_LMIN", "-8")), log_lambda_max=float(os.environ.get("HAPC_LMAX", "2")), grid_length=int(os.environ.get("HAPC_GRIDN", "12")))
MODEL_NAME = {"2": "PCHAR", "1": "PCHAL", "sv": "sv", "ridge": "ridge"}
rows = []

def design_file(stem, design):
    return D + stem + ("_raw.csv" if design == "raw" else ".csv")

def x_cols(df):
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

def record(**kw):
    kw["model"] = MODEL_NAME.get(kw["model"], kw["model"])
    rows.append(kw)

for design in ([] if os.environ.get("HAPC_SUMMARY_ONLY", "0") == "1" else DESIGNS):
    # ── (a) in-fill ─────────────────────────────────────────────────────────
    for target in ("level", "prev"):
        df = pd.read_csv(design_file(f"infill_ghana_child_iron_{target}", design))
        X = df[x_cols(df)].to_numpy(float); y = df["y_mod"].to_numpy(float); ynat = df["y_nat"].to_numpy(float)
        for norm in NORMS + ["ridge"]:
            rhos, lams, secs = [], [], []
            for rep in range(1, REPS + 1):
                folds = df[f"fold_rep{rep}"].to_numpy(int); pred = np.full(len(y), np.nan)
                for f in np.unique(folds):
                    te = folds == f; tr = ~te
                    if norm == "ridge":
                        pred[te] = ridge_predict(X[tr], y[tr], X[te]); lam, sec = np.nan, 0.0
                    else:
                        pred[te], lam, sec = fit_predict(X[tr], y[tr], X[te], norm)
                    lams.append(lam); secs.append(sec)
                rhos.append(rho(ynat, pred))
            record(design=design, estimand="infill", country="Ghana", outcome="child_iron", cell="Ghana child_iron", target=target, model=norm, max_degree=MAXDEG, n=len(y), p=X.shape[1],
                   spearman_median=float(np.median(rhos)), spearman_min=float(np.min(rhos)), spearman_max=float(np.max(rhos)),
                   median_lambda=float(np.nanmedian(lams)) if len(lams) else np.nan, seconds_per_fit=float(np.mean(secs)))
            print(f"[{design:3s}] infill {target:5s} {MODEL_NAME.get(norm, norm):6s} p={X.shape[1]:4d} rho median {np.median(rhos):+.3f} (draws {np.min(rhos):+.3f}..{np.max(rhos):+.3f}) | {np.mean(secs):.1f} s per fit", flush=True)

    # ── (b) leave-one-country-out ───────────────────────────────────────────
    for on in ("child_iron", "child_vitA"):
        for target in ("level", "prev"):
            df = pd.read_csv(design_file(f"loco_{on}_{target}", design))
            X = df[x_cols(df)].to_numpy(float); y = df["y_mod"].to_numpy(float); ynat = df["y_nat"].to_numpy(float); ctry = df["country"].to_numpy()
            for norm in NORMS + ["ridge"]:
                out = {}
                for cn in np.unique(ctry):
                    te = ctry == cn; tr = ~te
                    if norm == "ridge":
                        pred = ridge_predict(X[tr], y[tr], X[te]); lam, sec = np.nan, 0.0
                    else:
                        pred, lam, sec = fit_predict(X[tr], y[tr], X[te], norm)
                    out[cn] = rho(ynat[te], pred)
                    record(design=design, estimand="country", country=cn, outcome=on, cell=f"{cn} {on}", target=target, model=norm, max_degree=MAXDEG, n=int(te.sum()), p=X.shape[1],
                           spearman_median=float(out[cn]), spearman_min=np.nan, spearman_max=np.nan, median_lambda=float(lam), seconds_per_fit=float(sec))
                print(f"[{design:3s}] LOCO {on:10s} {target:5s} {MODEL_NAME.get(norm, norm):6s} p={X.shape[1]:4d} " + " | ".join(f"{k} {v:+.2f}" for k, v in out.items()) + f" | mean {np.mean(list(out.values())):+.3f}", flush=True)

if os.environ.get("HAPC_SUMMARY_ONLY", "0") == "1":
    res = pd.read_csv(D + OUT)
else:
    res = pd.DataFrame(rows)
    res.to_csv(D + OUT, index=False)

# ── paired summary: each model under each design, next to the domain index ──
# The index reference is matched on (country, outcome, target, estimand), so the
# in-fill row is Ghana's own cell and the LOCO row averages the held-out countries
try:
    cells = pd.read_csv(ROOT + "results/tables/protocol_v2/benchmarks_v2_cells.csv")
    idx = cells[cells["arm"] == "domain_index"][["country", "outcome", "target", "estimand", "spearman"]].rename(columns={"spearman": "domain_index"})
    res = res.merge(idx, on=["country", "outcome", "target", "estimand"], how="left")
except Exception as e:  # reference optional
    print("domain index reference unavailable:", e); res["domain_index"] = np.nan
summ = res.groupby(["estimand", "outcome", "target", "design", "model"], as_index=False)[["spearman_median", "domain_index"]].mean()
wide = summ.pivot_table(index=["estimand", "outcome", "target"], columns=["model", "design"], values="spearman_median")
wide.columns = [f"{m}_{d}" for m, d in wide.columns]
wide = wide.reset_index().merge(summ.groupby(["estimand", "outcome", "target"], as_index=False)["domain_index"].mean(), on=["estimand", "outcome", "target"])
wide.to_csv(D + OUT.replace(".csv", "_summary.csv"), index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
print(chr(10) + "== mean Spearman by cell (LOCO rows are means over the held-out countries) ==")
print(wide.round(3).to_string(index=False))
print("DONE", flush=True)
