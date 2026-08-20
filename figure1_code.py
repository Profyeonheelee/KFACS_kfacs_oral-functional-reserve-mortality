"""
figure1_code.py -- Analytic cohort flow counts for Figure 1B.
(Figure 1A is a hand-drawn conceptual schematic, not generated from data.)
Run 00_setup_and_mice.py first.
"""
import pandas as pd

pre = pd.read_pickle("final_preMI.pkl")
lm = pd.read_pickle("landmark_cohort.pkl")

n_total = len(pre)
n_natural_teeth = pre["natural_teeth"].notna().sum()
n_functional_dentition = pre["functional_dentition"].notna().sum()

mi = pd.read_pickle("mi_long_v2.pkl")
death_by_imp = mi.groupby("_imputation_id")["death_event"].sum()
followup_by_imp = mi.groupby("_imputation_id")["followup_years"].mean()

print("Figure 1B cohort-flow counts:")
print(f"  Analytic cohort (baseline dental records): N = {n_total:,}")
print(f"  Natural teeth available (pre-imputation):   n = {n_natural_teeth:,}")
print(f"  Functional dentition calculable:            n = {n_functional_dentition:,}")
print(f"  All-cause deaths (pooled across {len(death_by_imp)} imputations): "
      f"mean = {death_by_imp.mean():.0f} (range {death_by_imp.min():.0f}-{death_by_imp.max():.0f})")
print(f"  Mean follow-up: {followup_by_imp.mean():.1f} years")
print(f"  Wave 2 landmark cohort (alive + dental data at Wave 2): n = {len(lm):,}, "
      f"post-landmark deaths = {int(lm['post_landmark_event'].sum()):,}")
