"""
etable3_code.py -- Analytic cohort, variable availability, and death
ascertainment (pre-imputation completeness audit).
Run 00_setup_and_mice.py first.
"""
import pandas as pd

pre = pd.read_pickle("final_preMI.pkl")
N = len(pre)
def pct(n): return f"{n/N*100:.1f}"

avail = {
 "Natural teeth available": pre["natural_teeth"].notna().sum(),
 "Dental implant-supported restoration count available": pre["dental_implant_count"].notna().sum(),
 "Functional dentition calculable": pre["functional_dentition"].notna().sum(),
 "Mortality follow-up time available": pre["followup_years"].notna().sum(),
 "Age and sex available": pre[["age","sex"]].notna().all(axis=1).sum(),
 "Income level non-missing": ((pre["income2_1"].notna()) & (pre["income2_1"] != 9)).sum(),
 "Systemic disease count available": pre["systemic_disease_count"].notna().sum(),
 "CHS weight-loss component available": pre["chs_wtloss"].notna().sum(),
 "MNA screening risk available": pre["mna_scr_gr"].notna().sum(),
 "CHS frailty status calculable": pre["frailty_3cat"].notna().sum(),
 "MMSE-KC score available": pre["mmse_kc_score"].notna().sum(),
 "Systemic reserve adverse count available": pre["systemic_reserve_adverse_count"].notna().sum(),
}
print(f"Baseline dental records included in analytic dataset: {N:,} (100.0%)")
for label, n in avail.items():
    print(f"{label:55s} {n:,}  ({pct(n)}%)")

n_deaths = int(pre["death_event"].sum())
print(f"\nAll-cause death events: {n_deaths:,} ({pct(n_deaths)}%)")

dtm = pre["death_time_method"].value_counts()
print(f"\nDeath date quality:")
print(f"  Exact reported month: {dtm.get('reported_month',0):,} ({dtm.get('reported_month',0)/n_deaths*100:.1f}%)")
print(f"  Approximated (midpoint): {dtm.get('midpoint',0):,} ({dtm.get('midpoint',0)/n_deaths*100:.1f}%)")

print(
 "\nDeath ascertainment source (call/center/home visit) is unavailable in these source datasets: "
 "the per-wave response codes (center_pal_5 / home_pal_5 / call_pal_5) do not contain a death flag "
 "(code 4), so source-specific counts are not reported."
)
