"""
table2_code.py -- Cox proportional hazards models for all-cause mortality.
Model 1: age + sex + center + cohort start year
Model 2: Model 1 + education + income level + smoking + alcohol +
         systemic disease count + adverse systemic reserve count
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
    return d

imp_list = [build_covariates(mi[mi["_imputation_id"] == m]) for m in range(1, M + 1)]
CENTER_DUMMIES = [f"center_{i}" for i in range(2, 10)]
MODEL1_VARS = ["age", "female", "cohort2017"] + CENTER_DUMMIES
MODEL2_VARS = MODEL1_VARS + ["edu", "income_middle", "income_low", "smoking_history_ge100",
                              "ever_alcohol_use", "systemic_disease_count", "systemic_reserve_adverse_count"]

def fit_pool(exposure_term_fn, covar_vars):
    betas, variances = [], []
    for d in imp_list:
        dd = d.copy()
        dd["_exposure"] = exposure_term_fn(dd)
        cols = ["_exposure"] + covar_vars
        sub = dd[["followup_years", "death_event"] + cols].dropna()
        X = sub[cols].values.astype(float)
        res = fit_coxph(sub["followup_years"].values, sub["death_event"].values, X)
        betas.append(res["beta"]); variances.append(res["se"]**2)
    betas = np.array(betas); variances = np.array(variances)
    qbar, se, df = rubin_pool_vector(betas, variances)
    hr, lo, hi, p = pooled_summary(qbar, se, df)
    return dict(hr=hr[0], lo=lo[0], hi=hi[0], p=p[0], n=len(sub), n_events=int(sub["death_event"].sum()))

exposures = {
 "natural_teeth_per5":        lambda d: -d["natural_teeth"]/5,
 "natural_teeth_lt10":        lambda d: d["natural_teeth_lt10"],
 "functional_dentition_per5": lambda d: -d["functional_dentition"]/5,
 "functional_dentition_lt20": lambda d: d["functional_dentition_lt20"],
 "functional_dentition_lt10": lambda d: d["functional_dentition_lt10"],
 "implant_per2":              lambda d: d["dental_implant_count"]/2,
 "implant_any":               lambda d: d["any_dental_implant"],
}

def fmt_p(p): return "<.001" if p < 0.001 else f"{p:.3f}"[1:]
def fmt_hr(r): return f"{r['hr']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})"

rows = []
for name, fn in exposures.items():
    m1 = fit_pool(fn, MODEL1_VARS)
    m2 = fit_pool(fn, MODEL2_VARS)
    print(f"{name:28s} Model1 {fmt_hr(m1):>20s} P{fmt_p(m1['p'])}  n={m1['n']}/{m1['n_events']}")
    print(f"{'':28s} Model2 {fmt_hr(m2):>20s} P{fmt_p(m2['p'])}  n={m2['n']}/{m2['n_events']}")
    rows.append(dict(exposure=name, model="Model 1", **m1))
    rows.append(dict(exposure=name, model="Model 2", **m2))

pd.DataFrame(rows).to_csv("table2_output.csv", index=False)
print("\nSaved table2_output.csv")
