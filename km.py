import numpy as np
def kaplan_meier(time, event):
    """Simple Kaplan-Meier estimator (no CI band)."""
    time = np.asarray(time); event = np.asarray(event)
    order = np.argsort(time)
    time, event = time[order], event[order]
    times = np.unique(time[event == 1])
    surv = 1.0
    out_t, out_s = [0.0], [1.0]
    for t in times:
        at_risk = np.sum(time >= t)
        d = np.sum((time == t) & (event == 1))
        surv *= (1 - d/at_risk)
        out_t.append(t); out_s.append(surv)
    return np.array(out_t), np.array(out_s)
