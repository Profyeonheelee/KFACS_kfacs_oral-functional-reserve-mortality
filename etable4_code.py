"""
etable4_code.py -- Panel A (continuous distributions), Panel B (categorical
group sizes), Panel C (reclassification matrix) for baseline oral exposures.
Run 00_setup_and_mice.py first.
"""
import pandas as pd, numpy as np

mi = pd.read_pickle("mi_long_v2.pkl")
d1 = mi[mi["_imputation_id"] == 1].copy()
N = len(d1)

print("=== Panel A: continuous distributions ===")
for var in ["natural_teeth","dental_implant_count","functional_dentition"]:
    x = d1[var]
    print(f"{var:24s} n={N}  mean={x.mean():.1f}  sd={x.std():.1f}  "
          f"median={x.median():.0f}  IQR=[{x.quantile(.25):.0f},{x.quantile(.75):.0f}]  "
          f"min={x.min():.0f}  max={x.max():.0f}")

def cat4(x):
    if x >= 20: return ">=20"
    if x >= 10: return "10-19"
    if x >= 1: return "1-9"
    return "0"
def impcat(x):
    if x == 0: return "0"
    if x <= 2: return "1-2"
    if x <= 4: return "3-4"
    return ">=5"
d1["nt_cat"] = d1["natural_teeth"].apply(cat4)
d1["fd_cat"] = d1["functional_dentition"].apply(cat4)
d1["imp_cat"] = d1["dental_implant_count"].apply(impcat)

print("\n=== Panel B: categorical group sizes ===")
for label, col in [("Natural teeth","nt_cat"), ("Implant count","imp_cat"), ("Functional dentition","fd_cat")]:
    print(f"\n{label}:")
    print((d1[col].value_counts() / N * 100).round(1).to_string())

print("\n=== Panel C: reclassification matrix (rows=natural teeth group, cols=functional dentition group) ===")
matrix = pd.crosstab(d1["nt_cat"], d1["fd_cat"], margins=True)
print(matrix.reindex(index=[">=20","10-19","1-9","0","All"], columns=[">=20","10-19","1-9","0","All"]))
