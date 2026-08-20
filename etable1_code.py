"""
etable1_code.py -- Components of the adverse systemic reserve marker count,
by functional dentition group.
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np

mi = pd.read_pickle("mi_long_v2.pkl")
M = int(mi["_imputation_id"].max())
imp_list = [mi[mi["_imputation_id"] == m] for m in range(1, M + 1)]
pre = pd.read_pickle("final_preMI.pkl")

n_total = len(imp_list[0])
n_ge20 = int((imp_list[0]["functional_dentition_lt20"] == 0.0).sum())
n_lt20 = int((imp_list[0]["functional_dentition_lt20"] == 1.0).sum())

def pool_prop(cond_fn, group=None):
    props = []
    for d in imp_list:
        x = d if group is None else d[d["functional_dentition_lt20"] == group]
        props.append(cond_fn(x).mean())
    return np.mean(props)

comp_thresholds = {
 "low_albumin": lambda d: pd.Series(d["b_alb"] < 3.8),
 "anemia_low_hemoglobin": lambda d: pd.Series(np.where(d["sex"] == 2, d["b_hb"] < 12.0, d["b_hb"] < 13.0)),
 "elevated_hscrp": lambda d: pd.Series(d["b_hscrp"] > 3.0),
 "vitamin_d_deficiency": lambda d: pd.Series(d["b_25vtd"] < 20),
}
raw_col = {"low_albumin":"b_alb","anemia_low_hemoglobin":"b_hb","elevated_hscrp":"b_hscrp","vitamin_d_deficiency":"b_25vtd"}

for name, fn in comp_thresholds.items():
    o = pool_prop(fn); g0 = pool_prop(fn, 0.0); g1 = pool_prop(fn, 1.0)
    n_missing_orig = pre[raw_col[name]].isna().sum()
    print(f"{name:24s} present overall {o*100:5.1f}%  >=20 {g0*100:5.1f}%  <20 {g1*100:5.1f}%  "
          f"(n missing pre-imputation: {n_missing_orig})")

def pool_cont(var, group=None):
    means = [(d[var].dropna() if group is None else d.loc[d["functional_dentition_lt20"]==group, var].dropna()).mean()
             for d in imp_list]
    sd_avg = np.mean([(d.loc[d["functional_dentition_lt20"]==group, var].std() if group is not None else d[var].std())
                       for d in imp_list])
    return np.mean(means), sd_avg

o = pool_cont("systemic_reserve_adverse_count")
g0 = pool_cont("systemic_reserve_adverse_count", 0.0)
g1 = pool_cont("systemic_reserve_adverse_count", 1.0)
print(f"\nAdverse systemic reserve marker count: overall {o[0]:.1f}({o[1]:.1f})  "
      f">=20 {g0[0]:.1f}({g0[1]:.1f})  <20 {g1[0]:.1f}({g1[1]:.1f})")
print(f"\nN: total={n_total}  >=20={n_ge20}  <20={n_lt20}")
