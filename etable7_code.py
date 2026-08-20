"""
etable7_code.py -- Exploratory baseline implant-supported reclassification
analyses (natural teeth <20 reclassified to functional dentition >=20 via
implants, vs still <20, vs reference natural teeth >=20).
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np
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
    d["reclass_grp"] = np.select(
        [d["natural_teeth"] >= 20, (d["natural_teeth"] < 20) & (d["functional_dentition"] >= 20),
         (d["natural_teeth"] < 20) & (d["functional_dentition"] < 20)],
        ["ref", "reclass", "still_low"], default="none")
    d["any_dental_implant"] = (d["dental_implant_count"] >= 1).astype(float)
    return d

imp_list = [build_covariates(mi[mi["_imputation_id"] == m]) for m in range(1, M + 1)]
CENTER_DUMMIES = [f"center_{i}" for i in range(2, 10)]
MODEL1_VARS = ["age","female","cohort2017"] + CENTER_DUMMIES
MODEL2_VARS = MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
                              "systemic_disease_count","systemic_reserve_adverse_count"]

def fit_pool(exposure_fn, covars, filt=None):
    betas, variances = [], []
    for d in imp_list:
        dd = filt(d) if filt else d
        dd = dd.copy()
        dd["_exposure"] = exposure_fn(dd)
        cols = ["_exposure"] + covars
        sub = dd[["followup_years","death_event"] + cols].dropna()
        X = sub[cols].values.astype(float)
        res = fit_coxph(sub["followup_years"].values, sub["death_event"].values, X)
        betas.append(res["beta"]); variances.append(res["se"]**2)
    betas = np.array(betas); variances = np.array(variances)
    qbar, se, df = rubin_pool_vector(betas, variances)
    hr, lo, hi, p = pooled_summary(qbar, se, df)
    return dict(hr=hr[0], lo=lo[0], hi=hi[0], p=p[0], n=len(sub), n_events=int(sub["death_event"].sum()))

def fmt_p(p): return "<.001" if p < 0.001 else f"{p:.3f}"[1:]
def fmt_hr(r): return f"{r['hr']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})"

print("=== Reclassified / still-low vs reference (natural teeth >=20) ===")
for grp in ["reclass","still_low"]:
    filt = lambda d, grp=grp: d[d["reclass_grp"].isin(["ref", grp])]
    for mname, covars in [("M1",MODEL1_VARS),("M2",MODEL2_VARS)]:
        r = fit_pool(lambda d, grp=grp: (d["reclass_grp"]==grp).astype(float), covars, filt=filt)
        print(f"{grp:12s} {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

print("\n=== Within natural teeth <20 subgroup: reclassified vs still-low ===")
filt_lt20 = lambda d: d[d["reclass_grp"].isin(["reclass","still_low"])]
for mname, covars in [("M1",MODEL1_VARS),("M2",MODEL2_VARS)]:
    r = fit_pool(lambda d: (d["reclass_grp"]=="reclass").astype(float), covars, filt=filt_lt20)
    print(f"within_lt20 {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

print("\n=== Implant exposure within natural teeth <20 subgroup ===")
filt_nt20 = lambda d: d[d["natural_teeth"] < 20]
for mname, covars in [("M1",MODEL1_VARS),("M2",MODEL2_VARS)]:
    r = fit_pool(lambda d: d["any_dental_implant"], covars, filt=filt_nt20)
    print(f"any_implant {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")
for mname, covars in [("M1",MODEL1_VARS),("M2",MODEL2_VARS)]:
    r = fit_pool(lambda d: d["dental_implant_count"]/2, covars, filt=filt_nt20)
    print(f"implant_per2 {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")
