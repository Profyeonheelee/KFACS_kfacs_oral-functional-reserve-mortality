"""
Cox proportional hazards regression with Efron handling of tied event times,
implemented using Newton-Raphson optimization. Cumulative running sums are used
for nested risk sets ordered by descending event time.
"""
import numpy as np
from scipy import stats


def _prep(time, event, X):
    order = np.argsort(-time)  # descending time
    return time[order], event[order], X[order]


def _grad_hess_loglik(beta, time_desc, event_desc, X_desc):
    """time_desc/event_desc/X_desc must already be sorted by descending time."""
    n, p = X_desc.shape
    eta = X_desc @ beta
    w = np.exp(eta)

    grad = np.zeros(p)
    hess = np.zeros((p, p))
    loglik = 0.0

    S0 = 0.0
    S1 = np.zeros(p)
    S2 = np.zeros((p, p))

    i = 0
    while i < n:
        t = time_desc[i]
        j = i
        # gather the whole block of tied times (all rows with this exact time)
        while j < n and time_desc[j] == t:
            j += 1
        block = slice(i, j)
        # add this whole block to the risk-set running sums BEFORE using them,
        # since risk set for time t = all with time >= t (includes ties)
        w_block = w[block]
        X_block = X_desc[block]
        S0 += w_block.sum()
        S1 += (w_block[:, None] * X_block).sum(axis=0)
        S2 += (w_block[:, None, None] * X_block[:, :, None] * X_block[:, None, :]).sum(axis=0)

        death_mask = event_desc[block] == 1
        d = int(death_mask.sum())
        if d > 0:
            w_d = w_block[death_mask]
            X_d = X_block[death_mask]
            S0_d = w_d.sum()
            S1_d = (w_d[:, None] * X_d).sum(axis=0)
            S2_d = (w_d[:, None, None] * X_d[:, :, None] * X_d[:, None, :]).sum(axis=0)

            loglik += X_d.sum(axis=0) @ beta
            for l in range(d):
                frac = l / d
                S0_l = S0 - frac * S0_d
                S1_l = S1 - frac * S1_d
                S2_l = S2 - frac * S2_d
                Zbar = S1_l / S0_l
                loglik -= np.log(S0_l)
                grad += X_d[l] - Zbar
                hess += (S2_l / S0_l) - np.outer(Zbar, Zbar)
        i = j

    return grad, hess, loglik


def fit_coxph(time, event, X, max_iter=50, tol=1e-9):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    X = np.asarray(X, dtype=float)
    time_d, event_d, X_d = _prep(time, event, X)
    p = X.shape[1]
    beta = np.zeros(p)

    for _ in range(max_iter):
        grad, hess, _ = _grad_hess_loglik(beta, time_d, event_d, X_d)
        try:
            delta = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            delta = np.linalg.lstsq(hess, grad, rcond=None)[0]
        beta_new = beta + delta
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    grad, hess, loglik = _grad_hess_loglik(beta, time_d, event_d, X_d)
    cov = np.linalg.inv(hess)
    se = np.sqrt(np.diag(cov))
    z = beta / se
    pval = 2 * (1 - stats.norm.cdf(np.abs(z)))
    hr = np.exp(beta)
    ci_lo = np.exp(beta - 1.96 * se)
    ci_hi = np.exp(beta + 1.96 * se)

    return dict(beta=beta, se=se, cov=cov, hr=hr, ci_lo=ci_lo, ci_hi=ci_hi, z=z, p=pval,
                n=X.shape[0], n_events=int(event.sum()), loglik=loglik)
