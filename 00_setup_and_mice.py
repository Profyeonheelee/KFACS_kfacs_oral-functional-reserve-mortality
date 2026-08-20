"""
00_setup_and_mice.py
=====================
Builds the analytic dataset from the raw KFACS source files and runs
multiple imputation (m = 20). Run this once; every table/figure script
below loads its outputs (mi_long_v2.pkl, final_preMI.pkl, landmark_cohort.pkl,
physical_measures.pkl) instead of re-deriving them.

Cox regression is implemented in coxph.py using Efron handling of tied
event times and Newton-Raphson optimization. Multiple imputation is performed
with sklearn's IterativeImputer using BayesianRidge and posterior sampling.
"""
import pandas as pd, numpy as np, openpyxl
from pathlib import Path
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge

# ---- 0. Set this to the folder containing the raw source files ------------
DATA_DIR = Path("/path/to/KFACS_raw_data")
MAIN_FILE       = DATA_DIR / "_최종_2016_2024_KFACS.xlsx"
DENTAL_FILE     = DATA_DIR / "_최종_2016-2020_통합자료_치과_v20230904_LYH.xlsx"
TALKFILE        = DATA_DIR / "TalkFile_outcome_w1w5.xlsx"
FOLLOWUP_DETAIL = DATA_DIR / "KFACS_추적조사상세_최종_260608__2_.xlsx"

# =============================================================================
# 1. Baseline covariates (2016 & 2017 sheets of MAIN_FILE)
# =============================================================================
def extract_cols(path, sheet, targets, id_col="No"):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    wanted = set(targets) | {id_col}
    col_idx = {}
    for i, name in enumerate(header):
        if name in wanted and name not in col_idx:
            col_idx[name] = i
    data = []
    for row in rows:
        rid = row[col_idx[id_col]] if col_idx.get(id_col, -1) < len(row) else None
        if not rid:
            continue
        data.append({name: (row[idx] if idx < len(row) else None) for name, idx in col_idx.items()})
    wb.close()
    return pd.DataFrame(data)

BASELINE_TARGETS = [
 "center","sex","age",
 "dis1_1","dis1_2","dis1_3","dis1_4","dis1_5","dis1_6","dis1_7",
 "dis2_1","dis2_2","dis2_3","dis2_4",
 "dis3_1","dis3_2","dis3_3","dis3_4","dis3_5",
 "dis4_1","dis4_2","dis4_3",
 "dis5_1_1","dis5_1_2","dis5_2",
 "dis6_1","dis6_2",
 "dis7_1","dis7_2","dis7_3",
 "dis8_1","dis8_2","dis8_3","dis8_4","dis8_5","dis8_6","dis8_7","dis8_8","dis8_9","dis8_10_1","dis8_10_2","dis8_11",
 "dis9_1","dis9_2","dis9_3","dis9_4",
 "dis10_1","dis10_2_1","dis10_2_2",
 "smoke1","smoke2","smoke21_1","smoke22_1","smoke22_2","smoke23_1","smoke23_2","smoke23_3","smoke23_4","smoke23_5",
 "alcohl1","alcohl1_1","alcohl2_1","alcohl2_2",
 "edu","income1","income1_1","income1_2","income2","income2_1",
 "height","weight","chs_wtloss","mna_scr_gr",
 "b_alb","b_cre","b_hb","b_hscrp","b_25vtd",
 "chs1_1","chs1_2",
 "mmse_kc_score",
 "rt_mea1","lt_mea1","rt_mea2","lt_mea2","sppb1","sppb2","sppb","New_sppb_1","New_sppb_2",
 "panorama","panorama_op",
]

df16 = extract_cols(MAIN_FILE, "2016", BASELINE_TARGETS); df16["cohort_start_year"] = 2016
df17 = extract_cols(MAIN_FILE, "2017", BASELINE_TARGETS); df17["cohort_start_year"] = 2017
baseline = pd.concat([df16, df17], ignore_index=True).rename(columns={"No": "id"})
baseline.to_pickle("baseline.pkl")
print("baseline.pkl:", baseline.shape)

# ---- grip strength and gait speed variables used for descriptive summaries ----
for c in ["rt_mea1","lt_mea1","rt_mea2","lt_mea2"]:
    baseline[c] = pd.to_numeric(baseline[c], errors="coerce")
baseline["grip_strength_max"] = baseline[["rt_mea1","lt_mea1","rt_mea2","lt_mea2"]].max(axis=1)
avg_time = baseline[["New_sppb_1","New_sppb_2"]].mean(axis=1)
baseline["gait_speed_m_per_s"] = 4.0 / avg_time
baseline.loc[~baseline["gait_speed_m_per_s"].between(0.1, 3.0), "gait_speed_m_per_s"] = np.nan
baseline.loc[~baseline["grip_strength_max"].between(1, 80), "grip_strength_max"] = np.nan
baseline[["id","grip_strength_max","gait_speed_m_per_s"]].to_pickle("physical_measures.pkl")

# =============================================================================
# 2. Dental data (DENTAL_FILE, wide 3-block sheet: baseline / wave2 / wave3)
# =============================================================================
wb = openpyxl.load_workbook(DENTAL_FILE, read_only=True, data_only=True)
ws = wb["2016-2019 통합 panorama분석자료"]
rows = ws.iter_rows(values_only=True)
next(rows); next(rows)  # skip 2 header rows
BLOCKS = {"baseline_2016_2017": 0, "followup_2018_2019": 26, "wave2020": 52}
COLS = {"id":0,"center":1,"institution":2,"serial":3,"natural_teeth":4,"pontic":5,"root_rest":6,
        "missing_teeth":7,"decayed_teeth":8,"restored_teeth":9,"dental_implant_count":10,"wisdom_tooth":11,
        "sinus_mucosal_thickening":12,"pmi_ratio":13,"mi_distance":14,"gt_angle":15,"mci":16,"rabl":17,
        "etc":18,"panorama":19,"panorama_op":20}
dental_data = {k: [] for k in BLOCKS}
for row in rows:
    for block_name, offset in BLOCKS.items():
        rid = row[offset + COLS["id"]] if offset + COLS["id"] < len(row) else None
        if not rid:
            continue
        dental_data[block_name].append({name: (row[offset+idx] if offset+idx < len(row) else None)
                                         for name, idx in COLS.items()})
wb.close()
for name, recs in dental_data.items():
    pd.DataFrame(recs).to_pickle(f"dental_{name}.pkl")
    print(f"dental_{name}.pkl:", len(recs))

# =============================================================================
# 3. Derive baseline covariates, dental exposures, merge with curated outcome
# =============================================================================
dent_base = pd.read_pickle("dental_baseline_2016_2017.pkl")
d = dent_base[["id", "natural_teeth", "dental_implant_count"]].copy()
d["natural_teeth"] = pd.to_numeric(d["natural_teeth"], errors="coerce")
d["dental_implant_count"] = pd.to_numeric(d["dental_implant_count"], errors="coerce")
d.loc[~d["natural_teeth"].between(0, 32), "natural_teeth"] = np.nan
d.loc[d["dental_implant_count"] < 0, "dental_implant_count"] = np.nan
d["any_dental_implant"] = np.where(d["dental_implant_count"].isna(), np.nan, (d["dental_implant_count"] >= 1).astype(float))
d["functional_dentition"] = (d["natural_teeth"].fillna(0) + d["dental_implant_count"].fillna(0)).clip(upper=32)
d.loc[d["natural_teeth"].isna(), "functional_dentition"] = np.nan
d["functional_dentition_lt20"] = np.where(d["functional_dentition"].isna(), np.nan, (d["functional_dentition"] < 20).astype(float))
d["natural_teeth_lt10"] = np.where(d["natural_teeth"].isna(), np.nan, (d["natural_teeth"] < 10).astype(float))
d["functional_dentition_lt10"] = np.where(d["functional_dentition"].isna(), np.nan, (d["functional_dentition"] < 10).astype(float))

b = baseline.copy()
for labcol in ["b_alb","b_cre","b_hb","b_hscrp","b_25vtd","height","weight","edu","age","income2_1"]:
    b[labcol] = pd.to_numeric(b[labcol], errors="coerce")
b["low_albumin"] = np.where(b["b_alb"].isna(), np.nan, (b["b_alb"] < 3.8).astype(float))
b["anemia_low_hemoglobin"] = np.where(b["b_hb"].isna() | b["sex"].isna(), np.nan,
        np.where(b["sex"] == 2, b["b_hb"] < 12.0, b["b_hb"] < 13.0).astype(float))
b["elevated_hscrp"] = np.where(b["b_hscrp"].isna(), np.nan, (b["b_hscrp"] > 3.0).astype(float))
b["vitamin_d_deficiency"] = np.where(b["b_25vtd"].isna(), np.nan, (b["b_25vtd"] < 20).astype(float))
comp_cols = ["low_albumin","anemia_low_hemoglobin","elevated_hscrp","vitamin_d_deficiency"]
nonmiss = b[comp_cols].notna().sum(axis=1)
b["systemic_reserve_adverse_count"] = np.where(nonmiss >= 3, b[comp_cols].sum(axis=1, skipna=True), np.nan)
b["bmi"] = b["weight"] / (b["height"]/100.0)**2
b["smoking_history_ge100"] = b["smoke1"].map({1.0:0.0, 2.0:1.0, 3.0:0.0})
b["ever_alcohol_use"] = b["alcohl1"].map({1.0:0.0, 2.0:1.0})

# income_level categories based on income2_1 coding
# High: codes 1-2; Middle: codes 3-4; Low/no income: codes 5-6; code 9: missing
def income_cat_num(x):
    if pd.isna(x) or x == 9: return np.nan
    if x in (1,2): return 1   # High  (income2_1: >=300만원/month)
    if x in (3,4): return 2   # Middle (100-300만원)
    if x in (5,6): return 3   # Low/no income (<100만원 or none)
    return np.nan
b["income_level_num"] = b["income2_1"].apply(income_cat_num)

disease_items = ["dis1_1","dis1_2","dis1_3","dis1_4","dis1_5","dis1_6","dis1_7",
 "dis2_1","dis2_2","dis2_3","dis2_4","dis3_1","dis3_2","dis3_3","dis3_4","dis3_5",
 "dis4_1","dis4_2","dis4_3","dis5_1_1","dis5_1_2","dis5_2","dis6_1",
 "dis7_1","dis7_2","dis7_3","dis8_1","dis8_2","dis8_3","dis8_4","dis8_5","dis8_6",
 "dis8_7","dis8_8","dis8_9","dis8_10_1","dis8_11","dis9_1","dis9_2","dis9_3","dis9_4",
 "dis10_1","dis10_2_1"]  # dis6_2 (dementia) excluded from systemic disease count
present = [c for c in disease_items if c in b.columns]
disease_num = b[present].apply(pd.to_numeric, errors="coerce")
disease_num_bin = disease_num.where(disease_num.isin([0,1]))
b["systemic_disease_count"] = disease_num_bin.sum(axis=1, skipna=True)
b.loc[disease_num_bin.notna().sum(axis=1) == 0, "systemic_disease_count"] = np.nan

keep_cov = ["id","cohort_start_year","center","sex","age","edu","income2_1","income_level_num",
            "smoking_history_ge100","ever_alcohol_use","height","weight","bmi",
            "chs_wtloss","mna_scr_gr","b_alb","b_cre","b_hb","b_hscrp","b_25vtd",
            "systemic_reserve_adverse_count","systemic_disease_count","mmse_kc_score"]
b_final = b[keep_cov]

# ---- outcome and frailty variables imported from the wave-1 outcome file ----
talk = pd.read_excel(TALKFILE)
talk_bl = talk[talk["wave"] == 1].drop(columns=["wave"])

final = d.merge(b_final, on="id", how="left").merge(talk_bl, on="id", how="left")
final = final[final["id"] != "NO"]                     # exclude nonparticipant placeholder row
final = final.dropna(subset=["cohort_start_year"])      # require known cohort start year -> N = 2,715
final.to_pickle("final_preMI.pkl")
print("final_preMI.pkl:", final.shape)

# =============================================================================
# 4. Multiple imputation (m = 20)
#    Includes the primary exposure (natural_teeth/implant count), income level,
#    and death_event/followup_years. Missing values are handled within the
#    multiple-imputation workflow.
# =============================================================================
df = final.copy()
df.loc[df["chs_wtloss"] == 9, "chs_wtloss"] = np.nan
df.loc[df["income2_1"] == 9, "income2_1"] = np.nan

CONTINUOUS = ["age","edu","height","weight","bmi","b_alb","b_cre","b_hb","b_hscrp","b_25vtd",
              "systemic_reserve_adverse_count","systemic_disease_count","mmse_kc_score",
              "chs_total","adl_score","iadl_score","natural_teeth","dental_implant_count","followup_years"]
CATEGORICAL = ["sex","smoking_history_ge100","ever_alcohol_use","chs_wtloss","mna_scr_gr",
               "frailty_3cat","frailty_2cat","adl_2cat","adl_3cat","iadl_2cat","iadl_3cat",
               "income2_1","death_event"]
IMPUTE_VARS = CONTINUOUS + CATEGORICAL
VALID_RANGE = {"sex":(1,2),"smoking_history_ge100":(0,1),"ever_alcohol_use":(0,1),"chs_wtloss":(1,2),
 "mna_scr_gr":(1,3),"frailty_3cat":(0,2),"frailty_2cat":(0,1),"adl_2cat":(0,1),"adl_3cat":(0,2),
 "iadl_2cat":(0,1),"iadl_3cat":(0,2),"income2_1":(1,6),"death_event":(0,1)}
CLIP_CONT = {"natural_teeth":(0,32),"dental_implant_count":(0,30),"followup_years":(0.05,9.5),
 "age":(60,100),"bmi":(12,50),"b_alb":(1,6),"b_cre":(0.1,10),"b_hb":(5,20),"b_hscrp":(0,50),
 "b_25vtd":(0,100),"mmse_kc_score":(0,30),"edu":(0,25),"systemic_reserve_adverse_count":(0,4),
 "systemic_disease_count":(0,15),"chs_total":(0,5),"adl_score":(0,7),"iadl_score":(0,10)}

X = df[IMPUTE_VARS].copy()
M = 20
imputed_datasets = []
for m in range(1, M+1):
    imputer = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True, max_iter=10, random_state=m)
    Xi = pd.DataFrame(imputer.fit_transform(X), columns=IMPUTE_VARS, index=X.index)
    for c, (lo, hi) in VALID_RANGE.items():
        Xi[c] = Xi[c].round().clip(lo, hi)
    for c, (lo, hi) in CLIP_CONT.items():
        Xi[c] = Xi[c].clip(lo, hi)
    out = df.copy()
    for c in IMPUTE_VARS:
        out[c] = Xi[c].values
    out["any_dental_implant"] = (out["dental_implant_count"] >= 1).astype(float)
    out["functional_dentition"] = (out["natural_teeth"] + out["dental_implant_count"]).clip(upper=32)
    out["functional_dentition_lt20"] = (out["functional_dentition"] < 20).astype(float)
    out["natural_teeth_lt10"] = (out["natural_teeth"] < 10).astype(float)
    out["functional_dentition_lt10"] = (out["functional_dentition"] < 10).astype(float)
    out["income_level"] = out["income2_1"].apply(lambda x: "High" if x in (1,2) else ("Middle" if x in (3,4) else "Low/no income"))
    out["_imputation_id"] = m
    imputed_datasets.append(out)
    print(f"imputation {m}/{M} done")

mi_long = pd.concat(imputed_datasets, ignore_index=True)
mi_long.to_pickle("mi_long_v2.pkl")
print("mi_long_v2.pkl:", mi_long.shape)

# =============================================================================
# 5. Wave-2 landmark cohort (for eTable 5)
# =============================================================================
bl_age = pd.read_excel(FOLLOWUP_DETAIL, sheet_name="BL").rename(columns={"No.":"id"})
w2_age = pd.read_excel(FOLLOWUP_DETAIL, sheet_name="W2(2018-2019)").rename(columns={"No.":"id"})
age_merge = bl_age.merge(w2_age, on="id", how="inner")
age_merge["landmark_time"] = age_merge["W2_age_365"] - age_merge["BL_age_365"]

dent_fu = pd.read_pickle("dental_followup_2018_2019.pkl")
dent_fu["natural_teeth_w2"] = pd.to_numeric(dent_fu["natural_teeth"], errors="coerce")
dent_fu["implant_w2"] = pd.to_numeric(dent_fu["dental_implant_count"], errors="coerce")
dent_fu.loc[~dent_fu["natural_teeth_w2"].between(0,32), "natural_teeth_w2"] = np.nan
dent_fu.loc[dent_fu["implant_w2"] < 0, "implant_w2"] = np.nan
dent_fu["fd_w2"] = (dent_fu["natural_teeth_w2"].fillna(0) + dent_fu["implant_w2"].fillna(0)).clip(upper=32)
dent_fu.loc[dent_fu["natural_teeth_w2"].isna(), "fd_w2"] = np.nan
dent_fu = dent_fu[["id","natural_teeth_w2","implant_w2","fd_w2"]]

d1 = mi_long[mi_long["_imputation_id"] == 1].copy()  # representative imputation
lm = d1.merge(age_merge[["id","landmark_time"]], on="id", how="inner").merge(dent_fu, on="id", how="inner")
lm = lm.dropna(subset=["landmark_time","natural_teeth_w2","fd_w2"])
lm = lm[lm["followup_years"] >= lm["landmark_time"]].copy()
lm["post_landmark_followup"] = lm["followup_years"] - lm["landmark_time"]
lm["post_landmark_event"] = lm["death_event"]
lm["implant_increase"] = (lm["implant_w2"] > lm["dental_implant_count"]).astype(int)
lm["fd_increase"] = (lm["fd_w2"] > lm["functional_dentition"]).astype(int)
lm["fd_decrease"] = (lm["fd_w2"] < lm["functional_dentition"]).astype(int)
lm["new_implant_and_fd_gain"] = ((lm["implant_increase"]==1) & (lm["fd_increase"]==1)).astype(int)
lm["restored_to_fd20"] = ((lm["functional_dentition"] < 20) & (lm["fd_w2"] >= 20)).astype(int)
lm["baseline_fd_lt20"] = (lm["functional_dentition"] < 20).astype(int)
lm.to_pickle("landmark_cohort.pkl")
print("landmark_cohort.pkl:", lm.shape)

print("\nAll setup artifacts written. Proceed to the per-table/figure scripts.")
