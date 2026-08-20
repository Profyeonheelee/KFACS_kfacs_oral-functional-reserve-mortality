"""
efigure1_code.py -- Distribution of baseline oral functional reserve measures.
Run 00_setup_and_mice.py first.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mi = pd.read_pickle("mi_long_v2.pkl")
d1 = mi[mi["_imputation_id"] == 1]
N = len(d1)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
axes[0].hist(d1["natural_teeth"], bins=33, color="#4C72B0", edgecolor="white")
axes[0].set_title("Natural teeth"); axes[0].set_xlabel("Count"); axes[0].set_ylabel("Participants")
axes[1].hist(d1["dental_implant_count"], bins=28, color="#55A868", edgecolor="white")
axes[1].set_title("Implant-supported restorations"); axes[1].set_xlabel("Count")
axes[2].hist(d1["functional_dentition"], bins=33, color="#C44E52", edgecolor="white")
axes[2].set_title("Functional dentition"); axes[2].set_xlabel("Count")
plt.suptitle(f"eFigure 1. Distribution of Baseline Oral Functional Reserve Measures (N = {N:,})")
plt.tight_layout()
plt.savefig("efigure1_distributions.png", dpi=150)
print("Saved efigure1_distributions.png")
