"""
efigure3_code.py -- Sequential pathway attenuation of the low functional
dentition-mortality association, shown as violin (model-implied log-HR
distribution) + deterministic quantile scatter + trend line.
Run table3_code.py first (needs table3_output.csv).
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as st

t3 = pd.read_csv("table3_output.csv")
t3 = t3[t3["step"].str.match(r"^[1-5]\.|^6A\.")].reset_index(drop=True)  # steps 1-5 + 6A (fully adjusted)
labels = ["1. Demographic","2. +Socioeconomic/\nbehavioral","3. +Weight loss/\nnutrition","4. +Cognition",
          "5. +Systemic\ndisease","6A. +Frailty &\nreserve (full)"]

beta = np.log(t3["hr"].values)
se = (np.log(t3["hi"].values) - np.log(t3["lo"].values)) / (2*1.96)

fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(labels))
quantile_levels = [0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975]
rng_jitter = np.linspace(-0.12, 0.12, len(quantile_levels))

for i in range(len(labels)):
    grid = np.linspace(beta[i]-3.5*se[i], beta[i]+3.5*se[i], 200)
    dens = st.norm.pdf(grid, loc=beta[i], scale=se[i])
    dens = dens / dens.max() * 0.32
    hr_grid = np.exp(grid)
    ax.fill_betweenx(hr_grid, x[i]-dens, x[i]+dens, color="#4C72B0", alpha=0.25, lw=0)
    q_logHR = st.norm.ppf(quantile_levels, loc=beta[i], scale=se[i])
    ax.scatter(x[i]+rng_jitter, np.exp(q_logHR), color="#4C72B0", s=14, alpha=0.7, zorder=3)

ax.plot(x, t3["hr"], "-", color="black", lw=1.5, zorder=4)
ax.scatter(x, t3["hr"], color="black", s=45, zorder=5)
for i in range(len(labels)):
    p = 2*(1-st.norm.cdf(abs(beta[i]/se[i])))
    pfmt = "<.001" if p < 0.001 else f"P={p:.3f}"
    ax.annotate(f"HR={t3['hr'].iloc[i]:.2f}\n{pfmt}", (x[i], t3["hr"].iloc[i]),
                textcoords="offset points", xytext=(0,10), ha="center", fontsize=7.5, fontweight="bold")

ax.axhline(1.0, color="dimgray", ls="--", lw=1)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("Hazard ratio")
ax.set_title("eFigure 3. Sequential Pathway Attenuation of the\nLow Functional Dentition-Mortality Association")
ax.set_ylim(0.5, 2.6)
plt.tight_layout()
plt.savefig("efigure3_pathway.png", dpi=150)
print("Saved efigure3_pathway.png")
