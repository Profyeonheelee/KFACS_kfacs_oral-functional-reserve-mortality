"""
figure2_code.py -- Kaplan-Meier survival curves by natural teeth / functional
dentition categories. Uses imputation 1 of 20 as the representative dataset
(KM curves are not algebraically poolable across multiply-imputed datasets).
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from km import kaplan_meier   # see km.py in this package

mi = pd.read_pickle("mi_long_v2.pkl")
d1 = mi[mi["_imputation_id"] == 1].copy()

def cat4(x):
    if x >= 20: return ">=20"
    if x >= 10: return "10-19"
    if x >= 1: return "1-9"
    return "0"
d1["nt_cat"] = d1["natural_teeth"].apply(cat4)
d1["fd_cat"] = d1["functional_dentition"].apply(cat4)

lm = pd.read_pickle("landmark_cohort.pkl")  # for reclassification group sizes (Panel D)
d1["reclass_group"] = np.select(
    [d1["natural_teeth"] >= 20, (d1["natural_teeth"] < 20) & (d1["functional_dentition"] >= 20),
     (d1["natural_teeth"] < 20) & (d1["functional_dentition"] < 20)],
    ["Natural teeth >=20", "Reclassified (implants push to >=20)", "Natural teeth <20, still <20"], default="")

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel A: natural teeth categories
ax = axes[0,0]
for lv in [">=20","10-19","1-9","0"]:
    sub = d1[d1["nt_cat"] == lv]
    t, s = kaplan_meier(sub["followup_years"].values, sub["death_event"].values)
    ax.step(t, s, where="post", label=f"{lv} (n={len(sub)})")
ax.set_title("A. Natural teeth"); ax.set_xlabel("Years"); ax.set_ylabel("Survival probability"); ax.legend(fontsize=8)

# Panel B: functional dentition categories
ax = axes[0,1]
for lv in [">=20","10-19","1-9","0"]:
    sub = d1[d1["fd_cat"] == lv]
    t, s = kaplan_meier(sub["followup_years"].values, sub["death_event"].values)
    ax.step(t, s, where="post", label=f"{lv} (n={len(sub)})")
ax.set_title("B. Functional dentition"); ax.set_xlabel("Years"); ax.legend(fontsize=8)

# Panel C: severe tooth loss
ax = axes[1,0]
for lv, lab in [(True,"<10 teeth"), (False,">=10 teeth")]:
    sub = d1[(d1["natural_teeth"] < 10) == lv]
    t, s = kaplan_meier(sub["followup_years"].values, sub["death_event"].values)
    ax.step(t, s, where="post", label=f"{lab} (n={len(sub)})")
ax.set_title("C. Severe natural tooth loss"); ax.set_xlabel("Years"); ax.set_ylabel("Survival probability"); ax.legend(fontsize=8)

# Panel D: reclassification
ax = axes[1,1]
for lv in ["Natural teeth >=20","Reclassified (implants push to >=20)","Natural teeth <20, still <20"]:
    sub = d1[d1["reclass_group"] == lv]
    t, s = kaplan_meier(sub["followup_years"].values, sub["death_event"].values)
    ax.step(t, s, where="post", label=f"{lv} (n={len(sub)})")
ax.set_title("D. Functional reclassification"); ax.set_xlabel("Years"); ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("figure2_km.png", dpi=150)
print("Saved figure2_km.png")
