import numpy as np
from scipy import stats
def rubin_pool_vector(betas, variances):
    """betas: (m,p) array of point estimates; variances: (m,p) array of corresponding variances.
    Returns qbar (p,), se (p,), df (p,) via Rubin's rules (with Barnard-Rubin small-sample df)."""
    m = betas.shape[0]
    qbar = betas.mean(axis=0)
    ubar = variances.mean(axis=0)
    b = betas.var(axis=0, ddof=1)
    T = ubar + (1 + 1/m) * b
    se = np.sqrt(T)
    # Barnard-Rubin adjusted df
    lam = (1 + 1/m) * b / T
    lam = np.clip(lam, 1e-8, 1 - 1e-8)
    df_old = (m - 1) / lam**2
    return qbar, se, df_old
def pooled_summary(qbar, se, df, exponentiate=True):
    tcrit = stats.t.ppf(0.975, df)
    lo = qbar - tcrit * se
    hi = qbar + tcrit * se
    tstat = qbar / se
    p = 2 * (1 - stats.t.cdf(np.abs(tstat), df))
    if exponentiate:
        return np.exp(qbar), np.exp(lo), np.exp(hi), p
    return qbar, lo, hi, p
