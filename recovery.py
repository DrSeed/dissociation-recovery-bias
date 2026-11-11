# Tissue dissociation as a selective filter: who survives isolation, and how to see the bias.
import numpy as np
from scipy import stats


def survival_probability(fragility, gentleness, harshness):
    """Chance a cell survives dissociation and capture: falls with fragility, more steeply for harsh protocols."""
    return 1 / (1 + np.exp(-(gentleness - harshness * fragility)))


def fraction_with_ci(is_state, mouse):
    """State fraction with mice (not cells) as the replicates: mean of per-mouse fractions, t-based 95% CI."""
    per_mouse = np.array([is_state[mouse == m].mean() for m in np.unique(mouse)])
    mean = per_mouse.mean()
    half = stats.t.ppf(0.975, len(per_mouse) - 1) * per_mouse.std(ddof=1) / np.sqrt(len(per_mouse))
    return mean, mean - half, mean + half, per_mouse


def signature_score(X, genes):
    return X[:, genes].mean(axis=1)


def log_fold_change(X, group, genes):
    return X[group][:, genes].mean() - X[~group][:, genes].mean()
