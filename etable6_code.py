"""
etable6_code.py -- Sensitivity and extended Cox analyses: four-category
exposure models, dementia-exclusion sensitivity, complete-case pathway
model, and month-level death-date sensitivity.
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np
from coxph import fit_coxph
from pooling import rubin_pool_vector, pooled_summary
from scipy import stats as st

mi = pd.read_pickle("mi_long_v2.pkl")
M = int(mi["_imputation_id"].max())
baseline_raw = pd.read_pickle("baseline.pkl")[["id","dis6_2"]]
baseline_raw["dementia"] = pd.to_numeric(baseline_raw["dis6_2"], errors="coerce")

def build_covariates(d):
    d = d.copy()
    d["female"] = (d["sex"] == 2).astype(float)
    d["cohort2017"] = (d["cohort_start_year"] == 2017).astype(float)
    d["income_middle"] = (d["income_level"] == "Middle").astype(float)
    d["income_low"] = (d["income_level"] == "Low/no income").astype(float)
    for c in sorted(d["center"].dropna().unique()):
        if c == 1: continue
        d[f"center_{int(c)}"] = (d["center"] == c).astype(float)
    d = d.merge(baseline_raw[["id","dementia"]], on="id", how="left")
    def cat4(x):
        if x >= 20: return "ge20"
        if x >= 10: return "10-19"
        if x >= 1: return "1-9"
        return "0"
    d["nt_cat"] = d["natural_teeth"].apply(cat4)
    d["fd_cat"] = d["functional_dentition"].apply(cat4)
    d["natural_teeth_lt20"] = (d["natural_teeth"] < 20).astype(float)
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

print("=== Four-category exposure models (Model 2, reference = >=20) ===")
for var, cats in [("nt_cat", ["10-19","1-9","0"]), ("fd_cat", ["10-19","1-9","0"])]:
    for cat in cats:
        filt = lambda d, var=var, cat=cat: d[d[var].isin(["ge20", cat])]
        r = fit_pool(lambda d, var=var, cat=cat: (d[var]==cat).astype(float), MODEL2_VARS, filt=filt)
        print(f"{var}={cat:8s} {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

r = fit_pool(lambda d: d["natural_teeth_lt20"], MODEL2_VARS)
print(f"natural_teeth_lt20       {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")
r = fit_pool(lambda d: d["functional_dentition_lt20"], MODEL2_VARS)
print(f"functional_dentition_lt20 {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

print("\n=== Dementia-exclusion sensitivity ===")
excl = lambda d: d[d["dementia"] != 1]
r = fit_pool(lambda d: d["functional_dentition_lt10"], MODEL2_VARS, filt=excl)
print(f"functional_dentition_lt10, dementia excluded  {fmt_hr(r):>15s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")
r = fit_pool(lambda d: d["natural_teeth_lt10"], MODEL2_VARS, filt=excl)
print(f"natural_teeth_lt10, dementia excluded          {fmt_hr(r):>15s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

print("\n=== Month-level death-date sensitivity (excludes approximated/midpoint deaths) ===")
for exp_col, label in [("functional_dentition_lt10","functional_dentition_lt10"), ("natural_teeth_lt10","natural_teeth_lt10")]:
    filt_month = lambda d: d[~((d["death_event"]==1) & (d["death_time_method"]=="midpoint"))]
    r = fit_pool(lambda d, c=exp_col: d[c], MODEL2_VARS, filt=filt_month)
    print(f"{label:30s} {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

print("\n=== Complete-case fully adjusted pathway model (no MICE, listwise deletion) ===")
pre = pd.read_pickle("final_preMI.pkl")
pre["female"] = (pre["sex"] == 2).astype(float)
pre["cohort2017"] = (pre["cohort_start_year"] == 2017).astype(float)
for c in sorted(pre["center"].dropna().unique()):
    if c == 1: continue
    pre[f"center_{int(c)}"] = (pre["center"] == c).astype(float)
pre["income_middle"] = (pre["income_level_num"] == 2).astype(float)
pre["income_low"] = (pre["income_level_num"] == 3).astype(float)
pre.loc[pre["income_level_num"].isna(), ["income_middle","income_low"]] = np.nan
pre["chs_wtloss_yn"] = (pre["chs_wtloss"] == 1).astype(float)
pre["mna_risk"] = (pre["mna_scr_gr"] == 2).astype(float)
pre["mna_malnourished"] = (pre["mna_scr_gr"] == 3).astype(float)
pre.loc[pre["mna_scr_gr"].isna(), ["mna_risk","mna_malnourished"]] = np.nan
pre["weight_loss_flag"] = (pre["chs_wtloss"] == 1).astype(float)
pre["physical_frailty_burden_no_weight_loss"] = (pre["chs_total"] - pre["weight_loss_flag"]).clip(lower=0)
covars_cc = MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
    "chs_wtloss_yn","bmi","mna_risk","mna_malnourished","mmse_kc_score","systemic_disease_count",
    "physical_frailty_burden_no_weight_loss","systemic_reserve_adverse_count"]
cols = ["functional_dentition_lt10"] + covars_cc
sub = pre[["followup_years","death_event"] + cols].dropna()
X = sub[cols].values.astype(float)
res = fit_coxph(sub["followup_years"].values, sub["death_event"].values, X)
hr = np.exp(res["beta"][0]); lo=np.exp(res["beta"][0]-1.96*res["se"][0]); hi=np.exp(res["beta"][0]+1.96*res["se"][0])
p = 2*(1-st.norm.cdf(abs(res["z"][0])))
print(f"complete-case fully adjusted: HR={hr:.2f} ({lo:.2f}-{hi:.2f}) P{fmt_p(p)} n={len(sub)}/{int(sub['death_event'].sum())}")
