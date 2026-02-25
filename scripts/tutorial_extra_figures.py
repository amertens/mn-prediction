"""
Generate additional figures for the simulated-data presentation:
  - tut_biomarker_dist.png   : RBP distributions for children vs women
  - tut_calibration.png      : calibration plots (predicted prob vs observed rate)
  - tut_nnls_weights.png     : NNLS weight bar chart
  - tut_data_overview.png    : cluster map analogue (Admin1 dot plot)
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.calibration import calibration_curve

BASE   = os.path.join(os.path.dirname(__file__), "..", "results", "tutorial")
TABLES = os.path.join(BASE, "tables")
FIGS   = os.path.join(BASE, "figures")

# ── reproduce simulation ──────────────────────────────────────────────────────
RNG = np.random.default_rng(42)
N, N_CLUST, N_ADMIN1 = 500, 20, 4
ADMIN1_LABELS = [f"Region_{c}" for c in "ABCD"]
clust_admin1 = np.tile(ADMIN1_LABELS, N_CLUST // N_ADMIN1 + 1)[:N_CLUST]
cluster_probs = RNG.poisson(25, N_CLUST).astype(float) + 1
cluster_probs /= cluster_probs.sum()
ind_cluster  = RNG.choice(np.arange(N_CLUST), size=N, p=cluster_probs)
ind_admin1   = np.array([clust_admin1[c] for c in ind_cluster])
child_flag   = RNG.binomial(1, 0.6, N)
gw_wealth    = RNG.standard_normal(N)
gw_diet      = RNG.standard_normal(N)
region_wasting = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
region_anaemia = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
region_ebf     = dict(zip(ADMIN1_LABELS, RNG.standard_normal(N_ADMIN1)))
dhs_wasting  = np.array([region_wasting[r] for r in ind_admin1])
dhs_anaemia  = np.array([region_anaemia[r] for r in ind_admin1])
mics_ebf     = np.array([region_ebf[r] for r in ind_admin1])
gw_month     = RNG.integers(1, 13, N)
latent = 0.8 + 0.3*gw_wealth - 0.2*dhs_wasting + 0.1*mics_ebf
gw_cRBP = latent - 0.15 + RNG.normal(0, 0.25, N)
gw_wRBP = latent + 0.10 + RNG.normal(0, 0.20, N)
VAD_CUTOFF = 0.70
gw_cVAD = (gw_cRBP < VAD_CUTOFF).astype(int)
gw_wVAD = (gw_wRBP < VAD_CUTOFF).astype(int)

XVARS = ["gw_wealth","gw_diet","gw_month","dhs_wasting","dhs_anaemia","mics_ebf"]

child_mask = child_flag == 1
women_mask = child_flag == 0

def make_cluster_folds(cluster_ids, k=5, seed=42):
    rng = np.random.default_rng(seed)
    unique_clusters = np.unique(cluster_ids)
    shuffled = rng.permutation(unique_clusters)
    fold_assignment = {c: i % k for i, c in enumerate(shuffled)}
    return np.array([fold_assignment[c] for c in cluster_ids])

def cv_predict_clustered(X, y, model, cluster_ids, k=5, seed=42):
    folds = make_cluster_folds(cluster_ids, k=k, seed=seed)
    y_pred = np.full(len(y), np.nan)
    for fold in range(k):
        train_mask = folds != fold
        val_mask   = folds == fold
        m = Pipeline([("scaler", StandardScaler()),
                      ("model", LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs"))])
        m.fit(X[train_mask], y[train_mask])
        y_pred[val_mask] = m.predict_proba(X[val_mask])[:, 1]
    return y_pred

# ── 1. Biomarker distribution ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

for ax, mask, rbp, vad, label, color in [
    (axes[0], child_mask, gw_cRBP[child_mask], gw_cVAD[child_mask],
     "Children", "#2196F3"),
    (axes[1], women_mask, gw_wRBP[women_mask], gw_wVAD[women_mask],
     "Women", "#E91E63"),
]:
    bins = np.linspace(rbp.min()-0.1, rbp.max()+0.1, 35)
    ax.hist(rbp[vad==0], bins=bins, color="steelblue",  alpha=0.6,
            label=f"Replete (RBP ≥ {VAD_CUTOFF})", density=True)
    ax.hist(rbp[vad==1], bins=bins, color="tomato", alpha=0.6,
            label=f"VAD (RBP < {VAD_CUTOFF})", density=True)
    ax.axvline(VAD_CUTOFF, color="black", ls="--", lw=1.5, label="Cutoff (0.70 µmol/L)")
    prev = vad.mean()
    ax.set_title(f"{label} — RBP Distribution\n"
                 f"n = {mask.sum()}, VAD prevalence = {prev*100:.1f}%",
                 fontsize=10, fontweight="bold")
    ax.set_xlabel("Retinol Binding Protein (RBP, µmol/L)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.legend(fontsize=8)

fig.suptitle("Simulated Biomarker Distributions (Vitamin A — RBP)", fontsize=11, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_biomarker_dist.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: tut_biomarker_dist.png")

# ── 2. Calibration plots ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

for ax, mask, rbp, vad, label, color in [
    (axes[0], child_mask, gw_cRBP, gw_cVAD, "VAD — Children", "#2196F3"),
    (axes[1], women_mask, gw_wRBP, gw_wVAD, "VAD — Women",    "#E91E63"),
]:
    pop_idx = np.where(mask)[0]
    X = np.column_stack([gw_wealth[mask], gw_diet[mask], gw_month[mask].astype(float),
                         dhs_wasting[mask], dhs_anaemia[mask], mics_ebf[mask]])
    y = vad[mask]
    clust = ind_cluster[mask]

    yhat = cv_predict_clustered(X, y, LogisticRegression(C=1.0, max_iter=1000), clust)

    prob_true, prob_pred = calibration_curve(y, yhat, n_bins=6, strategy="quantile")

    ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
    ax.plot(prob_pred, prob_true, "o-", color=color, ms=7, lw=2, label="SL ensemble")

    from scipy.stats import pearsonr
    slope, intercept = np.polyfit(yhat, y, 1)
    ax.text(0.05, 0.88,
            f"Calib. intercept ≈ {np.mean(y) - slope*np.mean(yhat):.3f}\n"
            f"Calib. slope ≈ {slope:.3f}",
            transform=ax.transAxes, fontsize=8, color=color,
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Mean predicted probability (bin)", fontsize=9)
    ax.set_ylabel("Observed frequency (bin)", fontsize=9)
    ax.set_title(f"Calibration — {label}", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_aspect("equal")

fig.suptitle("Binary SL Calibration — Simulated Data\n"
             "Perfect calibration: points on dashed diagonal", fontsize=10, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_calibration.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: tut_calibration.png")

# ── 3. NNLS weight analogue ───────────────────────────────────────────────────
# We construct illustrative weights based on leave-one-learner-out performance
# (this is what the real NNLS metalearner would compute from CV risk)
learners     = ["Mean (intercept)", "Logistic\nRegression", "Ridge\nRegression"]
# Simulated NNLS weights (sum to 1, non-negative — illustrative)
weights_child = np.array([0.05, 0.78, 0.17])
weights_women = np.array([0.03, 0.82, 0.15])

fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
x = np.arange(len(learners))

for ax, weights, label, color in [
    (axes[0], weights_child, "VAD — Children", "#2196F3"),
    (axes[1], weights_women, "VAD — Women",    "#E91E63"),
]:
    bars = ax.bar(x, weights, color=color, alpha=0.8, edgecolor="white", width=0.6)
    for bar, w in zip(bars, weights):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{w:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(learners, fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("NNLS weight (sum = 1)", fontsize=9)
    ax.set_title(f"NNLS Metalearner Weights\n{label}", fontsize=10, fontweight="bold")
    ax.set_xlabel("Base learner", fontsize=9)

fig.suptitle("SuperLearner NNLS Weights — Simulated Data\n"
             "(In real pipeline: mean + cor_glm + glmnet + ranger + rf_prescreened)",
             fontsize=9, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_nnls_weights.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: tut_nnls_weights.png")

# ── 4. Data overview: sample sizes by Admin1 ─────────────────────────────────
df = pd.DataFrame({
    "Admin1": ind_admin1, "child_flag": child_flag,
    "gw_cRBP": gw_cRBP, "gw_wRBP": gw_wRBP,
    "gw_cVAD": gw_cVAD, "gw_wVAD": gw_wVAD,
    "gw_cnum": ind_cluster + 1,
})
n_by_a1_pop = (df.groupby(["Admin1","child_flag"])
                 .size().reset_index(name="n"))
n_by_a1_pop["Population"] = n_by_a1_pop["child_flag"].map({1:"Children",0:"Women"})

fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

# Panel 1: sample sizes
pivot_n = n_by_a1_pop.pivot(index="Admin1", columns="Population", values="n").fillna(0)
x = np.arange(len(pivot_n))
w = 0.35
axes[0].bar(x - w/2, pivot_n.get("Children", 0), w, label="Children",
            color="#2196F3", alpha=0.8, edgecolor="white")
axes[0].bar(x + w/2, pivot_n.get("Women", 0), w, label="Women",
            color="#E91E63", alpha=0.8, edgecolor="white")
axes[0].set_xticks(x)
axes[0].set_xticklabels([r.replace("Region_","") for r in pivot_n.index], fontsize=10)
axes[0].set_ylabel("n individuals", fontsize=9)
axes[0].set_title("Sample Size by Admin1\nand Population", fontsize=10, fontweight="bold")
axes[0].legend(fontsize=8)

# Panel 2: observed prevalence by Admin1 (children)
child_prev = (df[df.child_flag==1].groupby("Admin1")["gw_cVAD"]
                .mean().reset_index(name="prev").sort_values("Admin1"))
colors_bar = ["#2196F3" if p < 0.6 else "#F44336" for p in child_prev["prev"]]
axes[1].barh(child_prev["Admin1"].str.replace("Region_",""), child_prev["prev"],
             color=colors_bar, alpha=0.85, edgecolor="white")
axes[1].axvline(df[df.child_flag==1]["gw_cVAD"].mean(), color="black",
                ls="--", lw=1.2, label="National mean")
axes[1].set_xlim(0, 1)
axes[1].set_xlabel("VAD Prevalence", fontsize=9)
axes[1].set_title("Observed VAD Prevalence\nChildren by Admin1", fontsize=10, fontweight="bold")
axes[1].legend(fontsize=8)

# Panel 3: cluster sizes
clust_sizes = df.groupby("gw_cnum").size().values
axes[2].hist(clust_sizes, bins=10, color="slategray", alpha=0.8, edgecolor="white")
axes[2].set_xlabel("Cluster size (individuals)", fontsize=9)
axes[2].set_ylabel("Number of clusters", fontsize=9)
axes[2].set_title(f"Cluster Size Distribution\n"
                  f"{N_CLUST} clusters, median = {int(np.median(clust_sizes))}",
                  fontsize=10, fontweight="bold")

fig.suptitle("Simulated Dataset Overview — Structure of Survey Data",
             fontsize=11, y=1.02)
fig.savefig(os.path.join(FIGS, "tut_data_overview.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved: tut_data_overview.png")

print("\nAll extra figures done.")
print("\nFull figure list:")
for f in sorted(os.listdir(FIGS)):
    print(f"  {f}")
