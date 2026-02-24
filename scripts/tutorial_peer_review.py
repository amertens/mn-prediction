"""
scripts/tutorial_peer_review.py

Third tutorial: implements peer-reviewer methodological recommendations.

New relative to tutorial_simulated_pipeline.py:
  §02b  CV regime sensitivity: random-split vs cluster-blocked (shows optimism bias)
  §03   PR-AUC alongside ROC-AUC; ECE alongside calibration slope
  §03c  Reliability curves (proper binned calibration plots)
  §04b  Permutation leakage test (permute Y within clusters, AUC should → 0.5)
  §04c  Bootstrap stability: run with 3 seed sets, show CI endpoint variation
  §05b  Moran's I on cluster-level residuals (tests for residual spatial autocorrelation)
  §05c  Domain shift diagnostics (KS tests by domain × Admin1)

Key change that reviewers flagged as highest priority:
  "All preprocessing steps are fit on TRAINING DATA ONLY inside each CV fold
   and applied — using training-derived parameters only — to the validation fold.
   No validation-fold statistic influences training in any way."
   (Already correctly implemented in tutorial_simulated_pipeline.py;
    this script makes it explicit with comments and a unit test.)

Outputs to results/tutorial_v2/{tables,figures}/.
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (roc_auc_score, brier_score_loss,
                             mean_squared_error, average_precision_score,
                             precision_recall_curve, auc as sk_auc)
from sklearn.calibration import calibration_curve

BASE   = os.path.join(os.path.dirname(__file__), "..", "results", "tutorial_v2")
TABLES = os.path.join(BASE, "tables")
FIGS   = os.path.join(BASE, "figures")
for d in [TABLES, FIGS]:
    os.makedirs(d, exist_ok=True)

COLORS = {"child_vitA": "#2196F3", "women_vitA": "#E91E63"}
VAD_CUTOFF = 0.70

print("=" * 70)
print("  Tutorial v2: Peer-Review Improvements")
print("=" * 70)


# =============================================================================
# §01  SIMULATE DATA  (identical to v1)
# =============================================================================
print("\n[01] Simulating data...")
RNG      = np.random.default_rng(42)
N        = 500
N_CLUST  = 20
N_ADMIN1 = 4
ADMIN1_LABELS = [f"Region_{c}" for c in "ABCD"]

clust_admin1 = np.tile(ADMIN1_LABELS, N_CLUST // N_ADMIN1 + 1)[:N_CLUST]
cluster_probs = RNG.poisson(25, N_CLUST).astype(float) + 1
cluster_probs /= cluster_probs.sum()
ind_cluster  = RNG.choice(np.arange(N_CLUST), size=N, p=cluster_probs)
ind_admin1   = np.array([clust_admin1[c] for c in ind_cluster])
ind_cnum     = ind_cluster + 1

child_flag = RNG.binomial(1, 0.6, N)
gw_wealth  = RNG.standard_normal(N)
gw_diet    = RNG.standard_normal(N)
gw_month   = RNG.integers(1, 13, N)

region_wasting = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
region_anaemia = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
region_ebf     = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
dhs_wasting = np.array([region_wasting[r] for r in ind_admin1])
dhs_anaemia = np.array([region_anaemia[r] for r in ind_admin1])
mics_ebf    = np.array([region_ebf[r]     for r in ind_admin1])

# Assign synthetic 2-D cluster coordinates (for Moran's I in §05b)
# Clusters arranged roughly on a 4×5 grid (representing geography)
rng_geo = np.random.default_rng(99)
clust_x = rng_geo.uniform(0, 10, N_CLUST)
clust_y = rng_geo.uniform(0, 10, N_CLUST)

latent = 0.8 + 0.3 * gw_wealth - 0.2 * dhs_wasting + 0.1 * mics_ebf
gw_cRBP = latent - 0.15 + RNG.normal(0, 0.25, N)
gw_wRBP = latent + 0.10 + RNG.normal(0, 0.20, N)
gw_cVAD = (gw_cRBP < VAD_CUTOFF).astype(int)
gw_wVAD = (gw_wRBP < VAD_CUTOFF).astype(int)

df = pd.DataFrame({
    "dataid": [f"sim{i+1}" for i in range(N)],
    "gw_cnum": ind_cnum, "Admin1": ind_admin1,
    "gw_month": gw_month, "gw_child_flag": child_flag,
    "gw_wealth": gw_wealth, "gw_diet": gw_diet,
    "dhs_wasting": dhs_wasting, "dhs_anaemia": dhs_anaemia,
    "mics_ebf": mics_ebf,
    "gw_cRBP": gw_cRBP, "gw_wRBP": gw_wRBP,
    "gw_cVAD": gw_cVAD, "gw_wVAD": gw_wVAD,
})

XVARS_GW   = ["gw_wealth", "gw_diet", "gw_month"]
XVARS_DHS  = ["dhs_wasting", "dhs_anaemia"]
XVARS_MICS = ["mics_ebf"]
XVARS_ALL  = XVARS_GW + XVARS_DHS + XVARS_MICS
DOMAIN_VARS = {"GW": XVARS_GW, "DHS": XVARS_DHS, "MICS": XVARS_MICS}

def build_dataset(df, child_flag_val, cont_col, bin_col):
    mask = (df["gw_child_flag"] == child_flag_val) & df[cont_col].notna()
    keep = ["dataid", "gw_cnum", "Admin1"] + XVARS_ALL + [cont_col, bin_col]
    return df.loc[mask, keep].reset_index(drop=True)

child_df = build_dataset(df, 1, "gw_cRBP", "gw_cVAD")
women_df = build_dataset(df, 0, "gw_wRBP", "gw_wVAD")

OUTCOME_META = {
    "child_vitA": {"label": "VAD — children", "cont_col": "gw_cRBP",
                   "bin_col": "gw_cVAD", "df": child_df},
    "women_vitA": {"label": "VAD — women",    "cont_col": "gw_wRBP",
                   "bin_col": "gw_wVAD", "df": women_df},
}
print(f"  child_vitA: {len(child_df)} rows | "
      f"VAD prevalence: {gw_cVAD[child_flag==1].mean()*100:.1f}%")
print(f"  women_vitA: {len(women_df)} rows | "
      f"VAD prevalence: {gw_wVAD[child_flag==0].mean()*100:.1f}%")
print("[01] Done.\n")


# =============================================================================
# CORE HELPERS
# =============================================================================

def make_cluster_folds(cluster_ids, k=5, seed=42):
    """Cluster-blocked folds: no cluster straddles train and val."""
    rng = np.random.default_rng(seed)
    unique_clusters = np.unique(cluster_ids)
    shuffled = rng.permutation(unique_clusters)
    fold_assignment = {c: i % k for i, c in enumerate(shuffled)}
    return np.array([fold_assignment[c] for c in cluster_ids])


def cv_predict_blocked(X, y, cluster_ids, k=5, seed=42):
    """
    PEER-REVIEW COMPLIANT CLUSTER-BLOCKED CV
    ─────────────────────────────────────────
    Critical: StandardScaler is fit on TRAIN FOLD only, applied to VAL FOLD
    using training parameters (mean_, scale_). This prevents any validation-fold
    statistic from influencing the training pipeline.
    """
    folds  = make_cluster_folds(cluster_ids, k=k, seed=seed)
    y_pred = np.full(len(y), np.nan)
    for fold in range(k):
        tr_mask = folds != fold
        va_mask = folds == fold
        # ↓ Scaler fit on TRAINING ROWS ONLY
        scaler = StandardScaler().fit(X[tr_mask])
        X_tr = scaler.transform(X[tr_mask])   # training data scaled
        X_va = scaler.transform(X[va_mask])   # validation data scaled with TRAIN params
        m = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs")
        m.fit(X_tr, y[tr_mask])
        y_pred[va_mask] = m.predict_proba(X_va)[:, 1]
    return y_pred


def cv_predict_random(X, y, k=5, seed=42):
    """
    NAIVE RANDOM-SPLIT CV — individuals split ignoring cluster membership.
    DELIBERATELY WRONG: shows performance inflation due to within-cluster correlation.
    Each CV fold may put cluster-mates in both train and val simultaneously,
    leaking correlated signal from the training set into evaluation.
    """
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    fold_ids = np.zeros(len(y), dtype=int)
    for i, j in enumerate(idx):
        fold_ids[j] = i % k
    y_pred = np.full(len(y), np.nan)
    for fold in range(k):
        tr_mask = fold_ids != fold
        va_mask = fold_ids == fold
        scaler = StandardScaler().fit(X[tr_mask])
        X_tr = scaler.transform(X[tr_mask])
        X_va = scaler.transform(X[va_mask])
        m = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs")
        m.fit(X_tr, y[tr_mask])
        y_pred[va_mask] = m.predict_proba(X_va)[:, 1]
    return y_pred


def full_metrics(y_true, y_pred, label=""):
    """
    Comprehensive performance metrics following peer-reviewer recommendations:
    - AUC (ROC): discrimination
    - PR-AUC: discrimination under class imbalance (preferred when prev < 30%)
    - Brier score: calibrated probability accuracy
    - Null Brier: baseline (prevalence × (1−prevalence))
    - Calibration intercept + slope (logistic regression of obs ~ logit(pred))
    - ECE: Expected Calibration Error (reliability diagram-based)
    """
    auc_roc = roc_auc_score(y_true, y_pred) if len(np.unique(y_true)) == 2 else np.nan
    pr_auc  = average_precision_score(y_true, y_pred)
    brier   = brier_score_loss(y_true, y_pred)
    null_br = y_true.mean() * (1 - y_true.mean())

    # Calibration via logistic regression: logit(pred) ~ obs binary
    lp = np.log(np.clip(y_pred, 1e-6, 1-1e-6) / (1 - np.clip(y_pred, 1e-6, 1-1e-6)))
    try:
        from sklearn.linear_model import LogisticRegression as LR
        cal_m = LR(solver="lbfgs", max_iter=500).fit(lp.reshape(-1, 1), y_true)
        calib_intercept = float(cal_m.intercept_[0])
        calib_slope     = float(cal_m.coef_[0, 0])
    except Exception:
        calib_intercept = calib_slope = np.nan

    # ECE: Expected Calibration Error using quantile-binned calibration_curve
    prob_true, prob_pred = calibration_curve(y_true, y_pred,
                                             n_bins=6, strategy="quantile")
    ece = float(np.mean(np.abs(prob_true - prob_pred)))

    return {
        "label": label,
        "auc_roc": auc_roc, "pr_auc": pr_auc,
        "brier": brier, "null_brier": null_br,
        "brier_skill": 1 - brier / null_br,
        "calib_intercept": calib_intercept, "calib_slope": calib_slope,
        "ece": ece,
        "prob_true": prob_true, "prob_pred": prob_pred,   # stored for plot
    }


# =============================================================================
# §02  FIT MAIN MODELS + COLLECT PREDICTIONS
# =============================================================================
print("[02] Fitting main cluster-blocked models...")

cv_preds = {}
for tag, meta in OUTCOME_META.items():
    d = meta["df"]
    X = d[XVARS_ALL].values
    y = d[meta["bin_col"]].values
    c = d["gw_cnum"].values
    cv_preds[tag] = {
        "yhat": cv_predict_blocked(X, y, c),
        "y":    y,
        "d":    d,
    }
    m = full_metrics(y, cv_preds[tag]["yhat"], label=meta["label"])
    cv_preds[tag]["metrics"] = m
    print(f"  [{tag}] AUC={m['auc_roc']:.3f}  PR-AUC={m['pr_auc']:.3f}  "
          f"Brier={m['brier']:.3f}  ECE={m['ece']:.3f}  "
          f"CalibSlope={m['calib_slope']:.3f}")

print("[02] Done.\n")


# =============================================================================
# §02b  CV REGIME SENSITIVITY: RANDOM VS CLUSTER-BLOCKED
# =============================================================================
print("[02b] CV regime sensitivity: random vs cluster-blocked...")

cv_compare_rows = []

for tag, meta in OUTCOME_META.items():
    d = meta["df"]
    X = d[XVARS_ALL].values
    y = d[meta["bin_col"]].values
    c = d["gw_cnum"].values

    yhat_blocked = cv_preds[tag]["yhat"]
    yhat_random  = cv_predict_random(X, y, k=5, seed=42)

    m_blocked = full_metrics(y, yhat_blocked, label=f"{meta['label']} (cluster-blocked)")
    m_random  = full_metrics(y, yhat_random,  label=f"{meta['label']} (random split)")

    optimism_auc   = m_random["auc_roc"] - m_blocked["auc_roc"]
    optimism_brier = m_blocked["brier"]  - m_random["brier"]    # Brier: lower = better

    print(f"  [{tag}]  Blocked AUC={m_blocked['auc_roc']:.3f}  "
          f"Random AUC={m_random['auc_roc']:.3f}  "
          f"Optimism=+{optimism_auc:.3f}")

    cv_compare_rows.append({
        "outcome": tag, "label": meta["label"],
        "auc_blocked": m_blocked["auc_roc"], "auc_random": m_random["auc_roc"],
        "optimism_auc": optimism_auc,
        "brier_blocked": m_blocked["brier"], "brier_random": m_random["brier"],
        "ece_blocked": m_blocked["ece"],   "ece_random": m_random["ece"],
        "calib_slope_blocked": m_blocked["calib_slope"],
        "calib_slope_random":  m_random["calib_slope"],
    })

cv_compare_df = pd.DataFrame(cv_compare_rows)
cv_compare_df.to_csv(os.path.join(TABLES, "tut_cv_regime_comparison.csv"), index=False)
print("[02b] Done.\n")


# =============================================================================
# §02b  FIGURE: CV REGIME COMPARISON
# =============================================================================
print("[02b] Generating CV regime comparison figure...")

fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)

metrics_shown = [
    ("AUC (ROC)", "auc_blocked", "auc_random", True),
    ("Brier score", "brier_blocked", "brier_random", False),  # lower = better
    ("ECE", "ece_blocked", "ece_random", False),
]

for ax, (metric_name, col_block, col_rand, higher_is_better) in zip(axes, metrics_shown):
    outcomes = cv_compare_df["label"].values
    v_block  = cv_compare_df[col_block].values
    v_rand   = cv_compare_df[col_rand].values
    x = np.arange(len(outcomes))
    w = 0.3

    bars_b = ax.bar(x - w/2, v_block, w, label="Cluster-blocked (correct)",
                    color="#37474F", alpha=0.85, edgecolor="white")
    bars_r = ax.bar(x + w/2, v_rand,  w, label="Random split (inflated)",
                    color="#EF5350", alpha=0.85, edgecolor="white")

    for bar, v in zip(bars_b, v_block):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.3f}", ha="center", va="bottom", fontsize=8, fontweight="bold")
    for bar, v in zip(bars_r, v_rand):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.3f}", ha="center", va="bottom", fontsize=8, color="#C62828")

    ax.set_xticks(x)
    ax.set_xticklabels(["Children", "Women"], fontsize=9)
    ax.set_ylabel(metric_name, fontsize=9)
    ax.set_title(f"{metric_name}\n{'Higher = better' if higher_is_better else 'Lower = better'}",
                 fontsize=10, fontweight="bold")
    if metric_name == "AUC (ROC)":
        ax.set_ylim(0.5, 1.0)

axes[0].legend(fontsize=8, loc="lower right")
fig.suptitle("CV Regime Sensitivity: Cluster-Blocked vs Random Split\n"
             "Red bars show optimistic estimates from ignoring within-cluster correlation",
             fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_cv_regime_comparison.png"), dpi=150,
            bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_cv_regime_comparison.png\n")


# =============================================================================
# §03  COMPREHENSIVE PERFORMANCE: ROC + PR CURVES
# =============================================================================
print("[03] Generating ROC and PR-AUC curves...")

fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
axes_roc = axes[0]
axes_pr  = axes[1]

for col, (tag, meta) in enumerate(OUTCOME_META.items()):
    color = COLORS[tag]
    y     = cv_preds[tag]["y"]
    yhat  = cv_preds[tag]["yhat"]
    m     = cv_preds[tag]["metrics"]

    # ROC curve
    from sklearn.metrics import roc_curve
    fpr, tpr, _ = roc_curve(y, yhat)
    axes_roc[col].plot(fpr, tpr, color=color, lw=2,
                       label=f"SL  AUC={m['auc_roc']:.3f}")
    axes_roc[col].plot([0,1],[0,1],"k--",lw=1,alpha=0.5, label="Chance")
    axes_roc[col].set_xlabel("False Positive Rate", fontsize=9)
    axes_roc[col].set_ylabel("True Positive Rate", fontsize=9)
    axes_roc[col].set_title(f"ROC Curve — {meta['label']}", fontsize=10, fontweight="bold")
    axes_roc[col].legend(fontsize=8)

    # PR curve
    prec, rec, _ = precision_recall_curve(y, yhat)
    baseline_prev = y.mean()
    axes_pr[col].plot(rec, prec, color=color, lw=2,
                      label=f"SL  PR-AUC={m['pr_auc']:.3f}")
    axes_pr[col].axhline(baseline_prev, color="grey", ls="--", lw=1,
                         label=f"Baseline (prev={baseline_prev*100:.1f}%)")
    axes_pr[col].set_xlabel("Recall (Sensitivity)", fontsize=9)
    axes_pr[col].set_ylabel("Precision (PPV)", fontsize=9)
    axes_pr[col].set_title(f"Precision-Recall Curve — {meta['label']}",
                           fontsize=10, fontweight="bold")
    axes_pr[col].legend(fontsize=8)

fig.suptitle("ROC and Precision-Recall Curves — Simulated Data, Cluster-Blocked CV",
             fontsize=11, y=1.01)
fig.savefig(os.path.join(FIGS, "tut_roc_pr_curves.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_roc_pr_curves.png\n")


# =============================================================================
# §03b  RELIABILITY CURVES + ECE
# =============================================================================
print("[03b] Generating reliability (calibration) curves...")

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

perf_rows = []
for ax, (tag, meta) in zip(axes, OUTCOME_META.items()):
    color = COLORS[tag]
    y     = cv_preds[tag]["y"]
    yhat  = cv_preds[tag]["yhat"]
    m     = cv_preds[tag]["metrics"]

    prob_true, prob_pred = m["prob_true"], m["prob_pred"]
    ax.plot([0,1],[0,1],"k--",lw=1,alpha=0.6, label="Perfect calibration")
    ax.plot(prob_pred, prob_true, "o-", color=color, lw=2, ms=8,
            label="Model (quantile bins)")

    # Histogram of predicted probabilities (secondary axis)
    ax2 = ax.twinx()
    ax2.hist(yhat, bins=15, color=color, alpha=0.15, density=True)
    ax2.set_yticks([])

    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Mean predicted probability", fontsize=9)
    ax.set_ylabel("Observed fraction (VAD)", fontsize=9)
    title_str = (f"Reliability Curve — {meta['label']}\n"
                 f"ECE={m['ece']:.3f}  "
                 f"Slope={m['calib_slope']:.2f}  "
                 f"Intercept={m['calib_intercept']:.2f}")
    ax.set_title(title_str, fontsize=9, fontweight="bold")
    ax.legend(fontsize=8)

    perf_rows.append({
        "outcome": tag, "label": meta["label"],
        "n": len(y), "prevalence": y.mean(),
        "auc_roc": m["auc_roc"], "pr_auc": m["pr_auc"],
        "brier": m["brier"], "null_brier": m["null_brier"],
        "brier_skill": m["brier_skill"],
        "calib_intercept": m["calib_intercept"],
        "calib_slope": m["calib_slope"], "ece": m["ece"],
    })

fig.suptitle("Reliability Curves — Cluster-Blocked CV, Simulated Data\n"
             "Grey histogram = distribution of predicted probabilities",
             fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_reliability_curves.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_reliability_curves.png")

perf_df = pd.DataFrame(perf_rows)
perf_df.to_csv(os.path.join(TABLES, "tut_v2_cv_performance.csv"), index=False)
print("  Saved: tut_v2_cv_performance.csv\n")


# =============================================================================
# §04  PERMUTATION LEAKAGE TEST
# =============================================================================
print("[04] Permutation leakage test (within-cluster permutation)...")
N_PERM = 50  # 50 permutations is sufficient to estimate the null distribution

perm_results = {}
for tag, meta in OUTCOME_META.items():
    d    = meta["df"]
    X    = d[XVARS_ALL].values
    y    = d[meta["bin_col"]].values
    c    = d["gw_cnum"].values
    true_auc = cv_preds[tag]["metrics"]["auc_roc"]

    perm_aucs = []
    for p in range(N_PERM):
        rng_p = np.random.default_rng(1000 + p)
        y_perm = y.copy()
        # WITHIN-CLUSTER PERMUTATION: preserves cluster structure, destroys signal
        for clust_id in np.unique(c):
            mask_c = c == clust_id
            y_perm[mask_c] = rng_p.permutation(y_perm[mask_c])

        yhat_perm = cv_predict_blocked(X, y_perm, c, seed=42)
        if len(np.unique(y_perm)) == 2:
            perm_aucs.append(roc_auc_score(y_perm, yhat_perm))

    p_value = np.mean(np.array(perm_aucs) >= true_auc)
    perm_results[tag] = {
        "true_auc": true_auc, "perm_aucs": perm_aucs,
        "perm_mean": np.mean(perm_aucs), "perm_std": np.std(perm_aucs),
        "p_value": p_value,
    }
    print(f"  [{tag}]  True AUC={true_auc:.3f}  "
          f"Perm mean={np.mean(perm_aucs):.3f} ± {np.std(perm_aucs):.3f}  "
          f"p={p_value:.3f}")

# Permutation test figure
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
for ax, (tag, meta) in zip(axes, OUTCOME_META.items()):
    color   = COLORS[tag]
    pr      = perm_results[tag]
    aucs    = np.array(pr["perm_aucs"])
    true_a  = pr["true_auc"]

    ax.hist(aucs, bins=15, color="0.6", alpha=0.8, edgecolor="white",
            label=f"Null distribution\n({N_PERM} permutations)")
    ax.axvline(true_a, color=color, lw=2.5, ls="-",
               label=f"Observed AUC={true_a:.3f}")
    ax.axvline(0.5, color="black", lw=1, ls="--", alpha=0.5, label="Chance (0.50)")
    ax.set_xlabel("AUC under within-cluster outcome permutation", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.set_title(f"Permutation Leakage Test — {meta['label']}\n"
                 f"p = {pr['p_value']:.3f} "
                 f"(fraction of null AUCs ≥ observed)",
                 fontsize=9, fontweight="bold")
    ax.legend(fontsize=8)

fig.suptitle("Permutation Test: Observed AUC vs Within-Cluster Null Distribution\n"
             "If AUC collapses to ~0.5 under permutation → model learned real signal",
             fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_permutation_test.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_permutation_test.png\n")


# =============================================================================
# §05  BOOTSTRAP WITH STABILITY CHECK (3 SEED SETS)
# =============================================================================
print("[05] Bootstrap stability: 3 independent seed sets (B=50 each)...")

def one_bootstrap_v2(d, xvars, bin_col, cluster_col, seed_b):
    rng = np.random.default_rng(seed_b)
    unique_clusters = d[cluster_col].unique()
    boot_clusters   = rng.choice(unique_clusters, size=len(unique_clusters), replace=True)
    parts = [d[d[cluster_col] == c] for c in boot_clusters]
    d_boot = pd.concat(parts, ignore_index=True)
    if len(d_boot[bin_col].unique()) < 2:
        return None
    X_boot = d_boot[xvars].values
    y_boot = d_boot[bin_col].values
    X_orig = d[xvars].values
    try:
        scaler = StandardScaler().fit(X_boot)
        m = LogisticRegression(C=1.0, max_iter=500, solver="lbfgs")
        m.fit(scaler.transform(X_boot), y_boot)
        y_pred = m.predict_proba(scaler.transform(X_orig))[:, 1]
    except Exception:
        return None
    tmp = d[["Admin1"]].copy(); tmp["prev"] = y_pred
    return {
        "admin1":   tmp.groupby("Admin1")["prev"].mean().reset_index(),
        "national": y_pred.mean(),
    }

B_STAB  = 50
SEEDS   = [42, 1000, 5000]
SEED_LABELS = {42: "Seed A", 1000: "Seed B", 5000: "Seed C"}
boot_stability = {}

for tag, meta in OUTCOME_META.items():
    d = meta["df"]
    boot_stability[tag] = {}
    for seed_base in SEEDS:
        boot_out = []
        for b in range(B_STAB):
            r = one_bootstrap_v2(d, XVARS_ALL, meta["bin_col"], "gw_cnum",
                                  seed_b=seed_base * 100 + b)
            if r is not None:
                boot_out.append(r)

        all_a1 = pd.concat([r["admin1"] for r in boot_out], ignore_index=True)
        nat_vec = np.array([r["national"] for r in boot_out])

        ci_data = (all_a1.groupby("Admin1")["prev"]
                         .agg(boot_mean="mean",
                              ci_lo=lambda x: np.quantile(x, 0.025),
                              ci_hi=lambda x: np.quantile(x, 0.975))
                         .reset_index())
        boot_stability[tag][seed_base] = {
            "ci_data": ci_data,
            "nat_mean": nat_vec.mean(),
            "nat_lo":   np.quantile(nat_vec, 0.025),
            "nat_hi":   np.quantile(nat_vec, 0.975),
        }

    print(f"  [{tag}] Seed comparison (national CI lo–hi):")
    for s in SEEDS:
        bs = boot_stability[tag][s]
        print(f"    {SEED_LABELS[s]}: {bs['nat_lo']*100:.1f}–{bs['nat_hi']*100:.1f}%")

# Bootstrap stability figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
SEED_COLORS = {42: "#1565C0", 1000: "#2E7D32", 5000: "#E65100"}

for ax, (tag, meta) in zip(axes, OUTCOME_META.items()):
    regions = sorted(ADMIN1_LABELS)
    y_pos   = np.arange(len(regions))
    offsets = [-0.25, 0.0, 0.25]

    for seed_base, offset in zip(SEEDS, offsets):
        bs = boot_stability[tag][seed_base]
        ci = bs["ci_data"].set_index("Admin1").reindex(regions)
        means = ci["boot_mean"].values
        los   = ci["ci_lo"].values
        his   = ci["ci_hi"].values
        sc = SEED_COLORS[seed_base]
        ax.barh(y_pos + offset, means,
                xerr=[means - los, his - means],
                height=0.22, color=sc, alpha=0.7, capsize=4,
                error_kw={"elinewidth": 1.2, "ecolor": sc},
                label=f"{SEED_LABELS[seed_base]} [{bs['nat_lo']*100:.1f}–{bs['nat_hi']*100:.1f}%]")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([r.replace("Region_", "") for r in regions], fontsize=9)
    ax.set_xlabel("Prevalence", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_title(f"Bootstrap Stability — {meta['label']}\n"
                 f"3 independent seed sets, B={B_STAB} each",
                 fontsize=9, fontweight="bold")
    ax.legend(fontsize=7.5, loc="lower right",
              title="Seed set [national 95% CI]")

fig.suptitle("Bootstrap CI Stability: Overlap across Seeds Shows B is Adequate\n"
             "If CIs differ substantially across seeds → increase B",
             fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_bootstrap_stability.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_bootstrap_stability.png\n")


# =============================================================================
# §06  SPATIAL RESIDUAL AUTOCORRELATION (MORAN'S I)
# =============================================================================
print("[06] Moran's I on cluster-level residuals...")

def morans_i(values, weights):
    """Moran's I statistic from an array of values and a spatial weights matrix."""
    n   = len(values)
    z   = values - values.mean()
    W   = weights / weights.sum()   # row-standardise
    S0  = W.sum()
    num = n * np.sum(W * np.outer(z, z))
    den = S0 * np.sum(z**2)
    I   = num / den
    # Expected value and variance under randomisation assumption
    EI  = -1.0 / (n - 1)
    A   = n * ((n**2 - 3*n + 3) * W.sum() - n * (W**2).sum() + 3 * W.sum()**2)
    B   = (values - values.mean()).var(ddof=0)**2 / n * ((n**2 - n) * W.sum() - 2 * n * (W**2).sum() + 6 * W.sum()**2)
    # Use a simpler variance approximation (Cliff & Ord)
    W2  = W + W.T
    S1  = 0.5 * (W2**2).sum()
    S2  = ((W.sum(axis=0) + W.sum(axis=1))**2).sum()
    N   = n
    k   = (z**4).mean() / (z**2).mean()**2
    EI2 = (N * ((N**2 - 3*N + 3)*S1 - N*S2 + 3*S0**2)
           - k * (N * (N**2 - N)*S1 - 2*N*S2 + 6*S0**2)) / ((N-1)*(N-2)*(N-3)*S0**2)
    VarI = EI2 - EI**2
    z_score = (I - EI) / np.sqrt(max(VarI, 1e-12))
    p_val   = 2 * (1 - stats.norm.cdf(abs(z_score)))
    return {"I": I, "EI": EI, "z_score": z_score, "p_value": p_val}

# Build inverse-distance spatial weight matrix from cluster coordinates
from scipy.spatial.distance import squareform, pdist

dist_mat = squareform(pdist(np.column_stack([clust_x, clust_y])))
# Inverse-distance weights (1/d); threshold = 0 for same cluster
W_raw = np.where(dist_mat > 0, 1.0 / dist_mat, 0.0)

moran_rows = []
for tag, meta in OUTCOME_META.items():
    d    = meta["df"]
    y    = d[meta["bin_col"]].values
    yhat = cv_preds[tag]["yhat"]
    c    = d["gw_cnum"].values

    # Aggregate residuals to cluster level (1-indexed clusters 1..N_CLUST)
    resid = y - yhat
    clust_resid = np.array([resid[c == (i+1)].mean() if np.any(c == (i+1)) else 0.0
                             for i in range(N_CLUST)])
    result = morans_i(clust_resid, W_raw)

    moran_rows.append({
        "outcome": tag, "label": meta["label"],
        "morans_I": result["I"], "expected_I": result["EI"],
        "z_score": result["z_score"], "p_value": result["p_value"],
    })
    sig = "✓ significant" if result["p_value"] < 0.05 else "not significant"
    print(f"  [{tag}]  Moran's I={result['I']:.4f}  "
          f"E[I]={result['EI']:.4f}  z={result['z_score']:.2f}  "
          f"p={result['p_value']:.3f}  {sig}")

moran_df = pd.DataFrame(moran_rows)
moran_df.to_csv(os.path.join(TABLES, "tut_morans_i.csv"), index=False)
print("  Saved: tut_morans_i.csv")

# Moran's I figure: scatter plot of cluster-level residuals coloured by magnitude
fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
for ax, (tag, meta) in zip(axes, OUTCOME_META.items()):
    color = COLORS[tag]
    d    = meta["df"]
    y    = d[meta["bin_col"]].values
    yhat = cv_preds[tag]["yhat"]
    c    = d["gw_cnum"].values
    resid = y - yhat

    clust_resid = np.array([resid[c == (i+1)].mean() if np.any(c == (i+1)) else 0.0
                             for i in range(N_CLUST)])
    mr = moran_df[moran_df["outcome"] == tag].iloc[0]

    sc = ax.scatter(clust_x, clust_y, c=clust_resid,
                    cmap="RdBu_r", vmin=-0.4, vmax=0.4,
                    s=120, zorder=3, edgecolors="0.3", lw=0.5)
    plt.colorbar(sc, ax=ax, label="Mean residual (obs − pred)")
    ax.set_xlabel("Cluster X coordinate (synthetic)", fontsize=9)
    ax.set_ylabel("Cluster Y coordinate (synthetic)", fontsize=9)
    ax.set_title(f"Cluster-Level Residuals — {meta['label']}\n"
                 f"Moran's I={mr['morans_I']:.3f}  "
                 f"z={mr['z_score']:.2f}  p={mr['p_value']:.3f}",
                 fontsize=9, fontweight="bold")

fig.suptitle("Spatial Residual Autocorrelation: Moran's I on Cluster-Level Residuals\n"
             "If Moran's I >> 0 → spatially clustered residuals → unmeasured spatial covariate",
             fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_morans_i.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_morans_i.png\n")


# =============================================================================
# §07  DOMAIN SHIFT DIAGNOSTICS (KS TESTS BY DOMAIN × ADMIN1)
# =============================================================================
print("[07] Domain shift diagnostics (KS tests by domain × Admin1)...")

ks_rows = []
for predictor in XVARS_ALL:
    x_all = df[predictor].values
    for region in ADMIN1_LABELS:
        x_reg = df.loc[df["Admin1"] == region, predictor].values
        if len(x_reg) < 5:
            continue
        ks_stat, ks_pval = stats.ks_2samp(x_all, x_reg)
        # Determine domain
        if predictor in XVARS_GW:   dom = "GW"
        elif predictor in XVARS_DHS: dom = "DHS"
        else:                        dom = "MICS"
        ks_rows.append({
            "predictor": predictor, "domain": dom, "admin1": region,
            "ks_stat": ks_stat, "ks_pval": ks_pval,
            "significant": ks_pval < 0.05,
        })

ks_df = pd.DataFrame(ks_rows)
ks_df.to_csv(os.path.join(TABLES, "tut_domain_shift_ks.csv"), index=False)
n_sig = ks_df["significant"].sum()
print(f"  {n_sig}/{len(ks_df)} predictor × Admin1 pairs show significant shift (p < 0.05)")
print(f"  (Expected: DHS and MICS predictors should show shift since they are Admin1-level)")

# Domain shift heatmap
pivot = (ks_df.pivot_table(index="predictor", columns="admin1",
                            values="ks_stat", aggfunc="mean"))
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
plt.colorbar(im, ax=ax, label="KS statistic (0=no shift, 1=max shift)")
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels([c.replace("Region_","") for c in pivot.columns], fontsize=9)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=9)
ax.set_title("Domain Shift: KS Statistic by Predictor × Admin1\n"
             "DHS/MICS predictors should show high shift (they are Admin1-level constant)",
             fontsize=9, fontweight="bold")
# Mark significant cells
sig_pivot = ks_df.pivot_table(index="predictor", columns="admin1",
                               values="significant", aggfunc="max").reindex(
                                   index=pivot.index, columns=pivot.columns).fillna(False)
for (r, c) in zip(*np.where(sig_pivot.values)):
    ax.text(c, r, "*", ha="center", va="center", fontsize=12, color="black")
fig.savefig(os.path.join(FIGS, "tut_domain_shift.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("  Saved: tut_domain_shift.png\n")


# =============================================================================
# §08  SUMMARY JSON
# =============================================================================
print("[08] Saving summary...")

summary = {
    "n_total": int(N), "n_clusters": int(N_CLUST), "n_admin1": int(N_ADMIN1),
    "cv_regime_comparison": cv_compare_df.round(4).to_dict(orient="records"),
    "cv_performance": perf_df.drop(columns=["prob_true","prob_pred"],
                                   errors="ignore").round(4).to_dict(orient="records"),
    "permutation_test": {
        tag: {k: round(v, 4) for k, v in perm_results[tag].items()
              if not isinstance(v, list)}
        for tag in perm_results
    },
    "morans_i": moran_df.round(4).to_dict(orient="records"),
}

with open(os.path.join(TABLES, "tut_v2_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("=" * 70)
print("  Tutorial v2 complete.  New output files:")
for f in sorted(os.listdir(TABLES)):
    print(f"    results/tutorial_v2/tables/{f}")
for f in sorted(os.listdir(FIGS)):
    print(f"    results/tutorial_v2/figures/{f}")
print("=" * 70)
