"""
table1_code.py -- Baseline characteristics by functional dentition group.
Run 00_setup_and_mice.py first (needs mi_long_v2.pkl, physical_measures.pkl
in the same working directory).
"""
import pandas as pd, numpy as np
from scipy import stats
from pooling import rubin_pool_vector

mi = pd.read_pickle("mi_long_v2.pkl")
phys = pd.read_pickle("physical_measures.pkl")
mi = mi.merge(phys, on="id", how="left")
M = int(mi["_imputation_id"].max())
imp_list = [mi[mi["_imputation_id"] == m] for m in range(1, M + 1)]

n_total = len(imp_list[0])
n_ge20 = int((imp_list[0]["functional_dentition_lt20"] == 0.0).sum())
n_lt20 = int((imp_list[0]["functional_dentition_lt20"] == 1.0).sum())

def pool_cont(var, group=None):
    means = []
    for d in imp_list:
        x = d[var].dropna() if group is None else d.loc[d["functional_dentition_lt20"] == group, var].dropna()
        means.append(x.mean())
    sd_avg = np.mean([(d.loc[d["functional_dentition_lt20"] == group, var].std() if group is not None else d[var].std())
                       for d in imp_list])
    return np.mean(means), sd_avg

def pool_prop(cond_fn, group=None):
    props = []
    for d in imp_list:
        x = d if group is None else d[d["functional_dentition_lt20"] == group]
        props.append(cond_fn(x).mean())
    return np.mean(props)

def diff_pool_continuous(var):
    ests, varis = [], []
    for d in imp_list:
        x1 = d.loc[d["functional_dentition_lt20"] == 0.0, var].dropna()
        x2 = d.loc[d["functional_dentition_lt20"] == 1.0, var].dropna()
        ests.append(x1.mean() - x2.mean())
        varis.append(x1.var()/len(x1) + x2.var()/len(x2))
    ests, varis = np.array(ests), np.array(varis)
    qbar, se, df = rubin_pool_vector(ests.reshape(-1,1), varis.reshape(-1,1))
    return 2*(1-stats.t.cdf(np.abs(qbar[0]/se[0]), df[0]))

def diff_pool_binary(cond_fn):
    ests, varis = [], []
    for d in imp_list:
        x1 = cond_fn(d[d["functional_dentition_lt20"] == 0.0])
        x2 = cond_fn(d[d["functional_dentition_lt20"] == 1.0])
        p1, p2 = x1.mean(), x2.mean()
        ests.append(p1 - p2)
        varis.append(p1*(1-p1)/len(x1) + p2*(1-p2)/len(x2))
    ests, varis = np.array(ests), np.array(varis)
    qbar, se, df = rubin_pool_vector(ests.reshape(-1,1), varis.reshape(-1,1))
    return 2*(1-stats.t.cdf(np.abs(qbar[0]/se[0]), df[0]))

def smd_cont(m1,s1,m2,s2): return (m1-m2)/np.sqrt((s1**2+s2**2)/2)
def smd_bin(p1,p2): return (p1-p2)/np.sqrt((p1*(1-p1)+p2*(1-p2))/2)
def fmt_p(p): return "<.001" if p < 0.001 else f"{p:.3f}"[1:]
def fmt(m,s): return f"{m:.1f} ({s:.1f})"

print(f"N total={n_total}  >=20 n={n_ge20}  <20 n={n_lt20}\n")

CONT_VARS = ["age","edu","bmi","mmse_kc_score","systemic_disease_count","systemic_reserve_adverse_count",
             "natural_teeth","dental_implant_count","functional_dentition","b_alb",
             "grip_strength_max","gait_speed_m_per_s","followup_years"]
rows = []
for v in CONT_VARS:
    o = pool_cont(v); g0 = pool_cont(v,0.0); g1 = pool_cont(v,1.0)
    smd = smd_cont(g0[0],g0[1],g1[0],g1[1]); p = diff_pool_continuous(v)
    rows.append(dict(variable=v, overall=fmt(*o), ge20=fmt(*g0), lt20=fmt(*g1), smd=round(smd,2), p=fmt_p(p)))
    print(f"{v:32s} overall {fmt(*o):>14s}  >=20 {fmt(*g0):>14s}  <20 {fmt(*g1):>14s}  SMD {smd:+.2f}  P{fmt_p(p)}")

BIN_SPECS = [("sex",2.0,"female"), ("smoking_history_ge100",1.0,"smoking"), ("ever_alcohol_use",1.0,"alcohol"),
             ("any_dental_implant",1.0,"any_implant"), ("chs_wtloss",1.0,"weight_loss"), ("death_event",1.0,"death")]
for var, lv, label in BIN_SPECS:
    o = pool_prop(lambda d,var=var,lv=lv: d[var]==lv)
    g0 = pool_prop(lambda d,var=var,lv=lv: d[var]==lv, 0.0)
    g1 = pool_prop(lambda d,var=var,lv=lv: d[var]==lv, 1.0)
    smd = smd_bin(g0,g1); p = diff_pool_binary(lambda d,var=var,lv=lv: d[var]==lv)
    print(f"{label:32s} overall {o*100:5.1f}%        >=20 {g0*100:5.1f}%        <20 {g1*100:5.1f}%        SMD {smd:+.2f}  P{fmt_p(p)}")

for cat, name in [(0.0,"MNA_normal"),(2.0,"MNA_risk"),(3.0,"MNA_malnourished")]:
    o = pool_prop(lambda d,cat=cat: d["mna_scr_gr"]==cat)
    print(f"mna_scr_gr={name:16s} overall {o*100:5.1f}%")

for cat, name in [(0.0,"Robust"),(1.0,"Prefrail"),(2.0,"Frail")]:
    o = pool_prop(lambda d,cat=cat: d["frailty_3cat"]==cat)
    g0 = pool_prop(lambda d,cat=cat: d["frailty_3cat"]==cat, 0.0)
    g1 = pool_prop(lambda d,cat=cat: d["frailty_3cat"]==cat, 1.0)
    print(f"CHS frailty={name:12s} overall {o*100:5.1f}%  >=20 {g0*100:5.1f}%  <20 {g1*100:5.1f}%")

for lv in ["High","Middle","Low/no income"]:
    o = pool_prop(lambda d,lv=lv: d["income_level"]==lv)
    g0 = pool_prop(lambda d,lv=lv: d["income_level"]==lv, 0.0)
    g1 = pool_prop(lambda d,lv=lv: d["income_level"]==lv, 1.0)
    print(f"income={lv:16s} overall {o*100:5.1f}% (n={round(o*n_total)})  >=20 {g0*100:5.1f}% (n={round(g0*n_ge20)})  <20 {g1*100:5.1f}% (n={round(g1*n_lt20)})")

pd.DataFrame(rows).to_csv("table1_output.csv", index=False)
print("\nSaved table1_output.csv")
