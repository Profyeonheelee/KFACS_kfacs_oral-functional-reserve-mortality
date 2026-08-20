"""
table3_code.py -- Sequential pathway attenuation models.
Exposure: functional dentition <10 vs >=10.
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
    d["weight_loss_flag"] = (d["chs_wtloss"] == 1).astype(float)
    # NOTE: chs_total (curated, wave-1) is the simple sum of 5 CHS components
    # (exhaustion, weakness, slowness, low activity, weight loss), so this
    # subtraction is an EXACT non-weight-loss component count, not an approximation.
    d["physical_frailty_burden_no_weight_loss"] = (d["chs_total"] - d["weight_loss_flag"]).clip(lower=0)
    return d

imp_list = [build_covariates(mi[mi["_imputation_id"] == m]) for m in range(1, M + 1)]
CENTER_DUMMIES = [f"center_{i}" for i in range(2, 10)]
MODEL1_VARS = ["age", "female", "cohort2017"] + CENTER_DUMMIES

def fit_pool_cox(exposure_col, covars):
    betas, variances = [], []
    for d in imp_list:
        cols = [exposure_col] + covars
        sub = d[["followup_years", "death_event"] + cols].dropna()
        X = sub[cols].values.astype(float)
        res = fit_coxph(sub["followup_years"].values, sub["death_event"].values, X)
        betas.append(res["beta"]); variances.append(res["se"]**2)
    betas = np.array(betas); variances = np.array(variances)
    qbar, se, df = rubin_pool_vector(betas, variances)
    hr, lo, hi, p = pooled_summary(qbar, se, df)
    return dict(hr=hr[0], lo=lo[0], hi=hi[0], p=p[0], n=len(sub), n_events=int(sub["death_event"].sum()))

steps = {
 "1. Demographic": MODEL1_VARS,
 "2. + Socioeconomic/behavioral": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use"],
 "3. + Weight loss/nutrition": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
                                                "chs_wtloss","bmi","mna_scr_gr"],
 "4. + Cognition": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
                                    "chs_wtloss","bmi","mna_scr_gr","mmse_kc_score"],
 "5. + Systemic disease": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
                                           "chs_wtloss","bmi","mna_scr_gr","mmse_kc_score","systemic_disease_count"],
 "6A. Frailty + reserve (fully adjusted)": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100",
      "ever_alcohol_use","chs_wtloss","bmi","mna_scr_gr","mmse_kc_score","systemic_disease_count",
      "physical_frailty_burden_no_weight_loss","systemic_reserve_adverse_count"],
 "6B. Reserve excluded (frailty only)": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100",
      "ever_alcohol_use","chs_wtloss","bmi","mna_scr_gr","mmse_kc_score","systemic_disease_count",
      "physical_frailty_burden_no_weight_loss"],
 "6C. Frailty excluded (reserve only)": MODEL1_VARS + ["edu","income_middle","income_low","smoking_history_ge100",
      "ever_alcohol_use","chs_wtloss","bmi","mna_scr_gr","mmse_kc_score","systemic_disease_count",
      "systemic_reserve_adverse_count"],
}

def fmt_p(p): return "<.001" if p < 0.001 else f"{p:.3f}"[1:]

rows = []
demo_logHR = None
for step_name, covars in steps.items():
    res = fit_pool_cox("functional_dentition_lt10", covars)
    logHR = np.log(res["hr"])
    if demo_logHR is None and step_name.startswith("1."):
        demo_logHR = logHR
    atten = np.nan if step_name.startswith("1.") else 100*(1 - logHR/demo_logHR)
    print(f"{step_name:42s} HR={res['hr']:.2f} ({res['lo']:.2f}-{res['hi']:.2f}) P{fmt_p(res['p'])} "
          f"atten={atten:.1f}%" if not np.isnan(atten) else f"{step_name:42s} HR={res['hr']:.2f} ({res['lo']:.2f}-{res['hi']:.2f}) P{fmt_p(res['p'])} (reference)")
    rows.append(dict(step=step_name, **res, pct_attenuation=atten))

pd.DataFrame(rows).to_csv("table3_output.csv", index=False)
print("\nSaved table3_output.csv")
