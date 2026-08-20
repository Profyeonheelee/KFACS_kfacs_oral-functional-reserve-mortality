"""
etable5_code.py -- Dental follow-up change analyses and landmark Cox models
(Wave 2, 2018-2019 landmark; Wave 3, 2020 later-follow-up subset).
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np
from coxph import fit_coxph
from scipy import stats as st

lm = pd.read_pickle("landmark_cohort.pkl")

def build_covariates(d):
    d = d.copy()
    d["female"] = (d["sex"] == 2).astype(float)
    d["cohort2017"] = (d["cohort_start_year"] == 2017).astype(float)
    d["income_middle"] = (d["income_level"] == "Middle").astype(float)
    d["income_low"] = (d["income_level"] == "Low/no income").astype(float)
    for c in sorted(d["center"].dropna().unique()):
        if c == 1: continue
        d[f"center_{int(c)}"] = (d["center"] == c).astype(float)
    d["natural_teeth_change"] = d["natural_teeth_w2"] - d["natural_teeth"]
    d["implant_gain"] = d["implant_w2"] - d["dental_implant_count"]
    return d

lm2 = build_covariates(lm)
CENTER_DUMMIES = [f"center_{i}" for i in range(2, 10)]
MODEL1 = ["age","female","cohort2017"] + CENTER_DUMMIES
MODEL2 = MODEL1 + ["edu","income_middle","income_low","smoking_history_ge100","ever_alcohol_use",
                    "systemic_disease_count","systemic_reserve_adverse_count","natural_teeth","dental_implant_count"]
MODEL3 = MODEL2 + ["natural_teeth_change"]
REDUCED_CLINICAL = MODEL1 + ["systemic_disease_count","systemic_reserve_adverse_count"]

def fit_landmark(d, exposure_col, covars):
    cols = [exposure_col] + covars
    sub = d[["post_landmark_followup","post_landmark_event"] + cols].dropna()
    X = sub[cols].values.astype(float)
    res = fit_coxph(sub["post_landmark_followup"].values, sub["post_landmark_event"].values, X)
    hr = np.exp(res["beta"][0]); lo = np.exp(res["beta"][0]-1.96*res["se"][0]); hi = np.exp(res["beta"][0]+1.96*res["se"][0])
    p = 2*(1-st.norm.cdf(abs(res["z"][0])))
    return dict(hr=hr, lo=lo, hi=hi, p=p, n=len(sub), n_events=int(sub["post_landmark_event"].sum()))

def fmt_p(p): return "<.001" if p < 0.001 else f"{p:.3f}"[1:]
def fmt_hr(r): return f"{r['hr']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})"

print(f"=== Panel A ===\nLandmark cohort (Wave 2, alive + dental data): n={len(lm2)}, "
      f"post-landmark deaths={int(lm2['post_landmark_event'].sum())}\n")

print("=== Panel B: change phenotype counts ===")
for phen, label in [("implant_increase","Implant increase"), ("fd_increase","Functional dentition increase"),
                     ("fd_decrease","Functional dentition decrease"), ("new_implant_and_fd_gain","New implant + FD gain")]:
    n = int(lm2[phen].sum()); N = len(lm2)
    dd = int(lm2.loc[lm2[phen]==1, "post_landmark_event"].sum())
    print(f"{label:35s} {n}/{N} ({n/N*100:.1f}%)  post-landmark deaths {dd}/{n}")
sub_lt20 = lm2[lm2["baseline_fd_lt20"] == 1]
n = int(sub_lt20["restored_to_fd20"].sum()); N = len(sub_lt20)
dd = int(sub_lt20.loc[sub_lt20["restored_to_fd20"]==1, "post_landmark_event"].sum())
print(f"{'Restored to FD>=20 (baseline FD<20)':35s} {n}/{N} ({n/N*100:.1f}%)  post-landmark deaths {dd}/{n}")

print("\n=== Panel C: landmark Cox models ===")
for exp_col, label in [("implant_increase","Implant increase"), ("fd_increase","FD increase"),
                        ("new_implant_and_fd_gain","New implant + FD gain"), ("restored_to_fd20","Restored to FD>=20")]:
    for mname, covars in [("M1",MODEL1),("M2",MODEL2)]:
        r = fit_landmark(lm2, exp_col, covars)
        print(f"{label:28s} {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")
r = fit_landmark(lm2, "implant_increase", MODEL3)
print(f"{'Implant increase':28s} M3: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

increasers = lm2[lm2["implant_increase"] == 1].copy()
increasers["implant_gain_per2"] = increasers["implant_gain"]/2
for mname, covars in [("M2",MODEL2),("M3",MODEL3)]:
    r = fit_landmark(increasers, "implant_gain_per2", covars)
    print(f"{'Implant gain per 2 (among increasers)':28s} {mname}: {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

r = fit_landmark(sub_lt20, "restored_to_fd20", REDUCED_CLINICAL)
print(f"{'Restored FD>=20, reduced clinical model':40s} {fmt_hr(r):>22s} P{fmt_p(r['p'])}  n={r['n']}/{r['n_events']}")

# Panel A row 3: later (2020, Wave 3) follow-up subset
bl_age = pd.read_excel("KFACS_추적조사상세_최종_260608__2_.xlsx", sheet_name="BL").rename(columns={"No.":"id"})
w3_age = pd.read_excel("KFACS_추적조사상세_최종_260608__2_.xlsx", sheet_name="W3(2020-2021)").rename(columns={"No.":"id"})
age_merge = bl_age.merge(w3_age, on="id", how="inner")
age_merge["landmark_time"] = age_merge["W3_age_365"] - age_merge["BL_age_365"]
dent2020 = pd.read_pickle("dental_wave2020.pkl")
dent2020["natural_teeth_w3"] = pd.to_numeric(dent2020["natural_teeth"], errors="coerce")
dent2020 = dent2020[dent2020["natural_teeth_w3"].notna()][["id"]]
mi = pd.read_pickle("mi_long_v2.pkl")
d1 = mi[mi["_imputation_id"] == 1]
lm3 = d1.merge(age_merge[["id","landmark_time"]], on="id", how="inner").merge(dent2020, on="id", how="inner")
lm3 = lm3.dropna(subset=["landmark_time"])
lm3 = lm3[lm3["followup_years"] >= lm3["landmark_time"]]
print(f"\nLater dental follow-up subset (2020, Wave 3): n={len(lm3)}, deaths={int(lm3['death_event'].sum())}")
