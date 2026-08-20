"""
efigure2_code.py -- Deficit-based dose-response association below the
functional dentition threshold (natural teeth AND functional dentition,
both expressed as deficit below 20 units).
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from coxph import fit_coxph
from pooling import rubin_pool_vector, pooled_summary

mi = pd.read_pickle("mi_long_v2.pkl")
M = int(mi["_imputation_id"].max())

def build_covariates(d):
    d = d.copy()
    d["female"] = (d["sex"] == 2).astype(float)
    d["cohort2017"] = (d["cohort_start_year"] == 2017).astype(float)
    d["income_middle"] = (d["income_level"] == "Middle").astype(float)
    d["income_low"] = (d["income_level"] == "Low/no income").astype(float)
    for c in sorted(d["center"].dropna().unique()):
        if c == 1: continue
        d[f"center_{int(c)}"] = (d["center"] == c).astype(float)
    d["nt_deficit_per5"] = np.maximum(20 - d["natural_teeth"], 0) / 5.0
    d["fd_deficit_per5"] = np.maximum(20 - d["functional_dentition"], 0) / 5.0
    return d

imp_list = [build_covariates(mi[mi["_imputation_id"] == m]) for m in range(1, M + 1)]
CENTER_DUMMIES = [f"center_{i}" for i in range(2, 10)]
MODEL2_VARS = ["age","female","cohort2017"] + CENTER_DUMMIES + ["edu","income_middle","income_low",
               "smoking_history_ge100","ever_alcohol_use","systemic_disease_count","systemic_reserve_adverse_count"]

results = {}
for expo in ["nt_deficit_per5", "fd_deficit_per5"]:
    betas, variances = [], []
    for d in imp_list:
        cols = [expo] + MODEL2_VARS
        sub = d[["followup_years","death_event"] + cols].dropna()
        X = sub[cols].values.astype(float)
        res = fit_coxph(sub["followup_years"].values, sub["death_event"].values, X)
        betas.append(res["beta"]); variances.append(res["se"]**2)
    betas = np.array(betas); variances = np.array(variances)
    qbar, se, df = rubin_pool_vector(betas, variances)
    hr, lo, hi, p = pooled_summary(qbar, se, df)
    results[expo] = dict(beta=qbar[0], se=se[0], hr=hr[0], lo=lo[0], hi=hi[0], p=p[0])
    print(f"{expo}: HR per 5-unit deficit = {hr[0]:.2f} (95% CI {lo[0]:.2f}-{hi[0]:.2f}), P={p[0]:.4f}")

deficits = np.arange(0, 21, 1)
fig, ax = plt.subplots(figsize=(7.5, 5.5))
colors = {"nt_deficit_per5": "#4C72B0", "fd_deficit_per5": "#C44E52"}
labels = {"nt_deficit_per5": "Natural teeth deficit", "fd_deficit_per5": "Functional dentition deficit"}
for key in ["nt_deficit_per5", "fd_deficit_per5"]:
    beta, se = results[key]["beta"], results[key]["se"]
    curve = np.exp(beta * deficits / 5)
    lo = np.exp((beta - 1.96*se) * deficits / 5)
    hi = np.exp((beta + 1.96*se) * deficits / 5)
    ax.plot(deficits, curve, color=colors[key], lw=2,
            label=f"{labels[key]} (HR per 5-unit = {results[key]['hr']:.2f}, "
                  f"95% CI {results[key]['lo']:.2f}-{results[key]['hi']:.2f})")
    ax.fill_between(deficits, lo, hi, color=colors[key], alpha=0.15)
ax.axhline(1.0, color="gray", ls="--", lw=1)
ax.set_xlabel("Deficit below 20 units (teeth or functional dentition)")
ax.set_ylabel("Adjusted hazard ratio for all-cause mortality")
ax.set_title("eFigure 2. Deficit-Based Dose-Response Association\nBelow the Functional Dentition Threshold")
ax.legend(fontsize=8, loc="upper left")
plt.tight_layout()
plt.savefig("efigure2_doseresponse.png", dpi=150)
print("Saved efigure2_doseresponse.png")
