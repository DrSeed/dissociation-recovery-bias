# Self-contained demo: injured cortex, two isolation protocols, four mice each; what the data says vs the tissue.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from recovery import survival_probability, fraction_with_ci, signature_score, log_fold_change

PROTOCOLS = {"gentle, low-shear": (2.5, 0.8, 0.1), "harsh, high-shear": (2.0, 3.0, 0.8)}  # (gentleness, harshness, stress)
GENES = {"homeostatic": range(0, 10), "reactive": range(10, 20), "remodelling": range(20, 30), "stress (IEG)": range(30, 38)}
N_GENES = 120


def tissue(n, reactive_share, rng):
    """True astrocytes in the tissue: reactive cells carry a remodelling programme of varying strength."""
    reactive = rng.uniform(size=n) < reactive_share
    remodel = np.where(reactive, rng.beta(2, 2, n), 0.0)
    fragility = rng.normal(0, 0.3, n) + np.where(reactive, 0.4 + 1.6 * remodel, 0)
    return reactive, remodel, fragility


def expression(reactive, remodel, stress, rng):
    X = rng.normal(0, 1, (len(reactive), N_GENES))
    X[:, list(GENES["homeostatic"])] += np.where(reactive, 0.3, 1.5)[:, None]
    X[:, list(GENES["reactive"])] += np.where(reactive, 1.5, 0.0)[:, None]
    X[:, list(GENES["remodelling"])] += (2.5 * remodel)[:, None]
    X[:, list(GENES["stress (IEG)"])] += stress[:, None]
    return X


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(12)
    true_share, n_tissue, n_mice = 0.30, 4000, 4

    # the tissue itself (what imaging or sc-qPCR of the intact tissue would tell you)
    t_reac, t_rem, _ = tissue(20000, true_share, rng)
    X_true = expression(t_reac, t_rem, np.zeros(len(t_reac)), rng)

    captured = {}
    for name, (gentle, harsh, stress_level) in PROTOCOLS.items():
        reac, rem, mouse, Xs = [], [], [], []
        for m in range(n_mice):
            r, mo, frag = tissue(n_tissue, true_share * rng.uniform(0.8, 1.2), rng)
            keep = rng.uniform(size=n_tissue) < survival_probability(frag, gentle, harsh)
            idx = rng.choice(np.flatnonzero(keep), 500, replace=False)
            stress = rng.gamma(2, stress_level / 2, len(idx))
            reac.append(r[idx]); rem.append(mo[idx]); mouse.append(np.full(len(idx), m))
            Xs.append(expression(r[idx], mo[idx], stress, rng))
        captured[name] = (np.concatenate(reac), np.concatenate(rem), np.concatenate(mouse), np.vstack(Xs))

    rows = [{"source": "tissue (truth)", "reactive_pct": 100 * t_reac.mean(), "ci_low": np.nan, "ci_high": np.nan,
             "remodelling_in_reactive": t_rem[t_reac].mean(),
             **{f"logFC_{g}": log_fold_change(X_true, t_reac, list(ix)) for g, ix in GENES.items()},
             "stress_score_mean": signature_score(X_true, list(GENES["stress (IEG)"])).mean()}]
    for name, (reac, rem, mouse, X) in captured.items():
        mean, lo, hi, _ = fraction_with_ci(reac, mouse)
        rows.append({"source": name, "reactive_pct": 100 * mean, "ci_low": 100 * lo, "ci_high": 100 * hi,
                     "remodelling_in_reactive": rem[reac].mean(),
                     **{f"logFC_{g}": log_fold_change(X, reac, list(ix)) for g, ix in GENES.items()},
                     "stress_score_mean": signature_score(X, list(GENES["stress (IEG)"])).mean()})
    table = pd.DataFrame(rows).round(3)
    table.to_csv("results/summary.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colours = {"gentle, low-shear": "#4c72b0", "harsh, high-shear": "#c44e52"}
    for k, (name, (reac, rem, mouse, X)) in enumerate(captured.items()):
        mean, lo, hi, per = fraction_with_ci(reac, mouse)
        axes[0].scatter(np.full(len(per), k) + rng.uniform(-0.08, 0.08, len(per)), 100 * per, color=colours[name])
        axes[0].errorbar(k + 0.2, 100 * mean, yerr=[[100 * (mean - lo)], [100 * (hi - mean)]], fmt="s", color="k", capsize=4)
    axes[0].axhline(100 * true_share, color="k", ls="--", label="true share in tissue")
    axes[0].set_xticks([0, 1], list(captured))
    axes[0].set_xlim(-0.5, 1.5)
    axes[0].set_ylabel("% reactive astrocytes (dots: mice; square: mean, 95% CI)")
    axes[0].set_title("Harsh isolation under-recovers reactive cells")
    axes[0].legend(fontsize=8)

    bins = np.linspace(0, 1, 26)
    axes[1].hist(t_rem[t_reac], bins=bins, density=True, histtype="step", lw=2, color="k", label="tissue (truth)")
    for name, (reac, rem, _, _) in captured.items():
        axes[1].hist(rem[reac], bins=bins, density=True, histtype="step", lw=2, color=colours[name], label=name)
    axes[1].set_xlabel("remodelling programme strength (reactive cells)")
    axes[1].set_ylabel("density")
    axes[1].set_title("...and loses the most remodelled ones first")
    axes[1].legend(fontsize=8)

    groups = list(GENES)
    x = np.arange(len(groups))
    for k, src in enumerate(["tissue (truth)"] + list(captured)):
        vals = [table.loc[table.source == src, f"logFC_{g}"].item() for g in groups]
        axes[2].bar(x + (k - 1) * 0.27, vals, 0.27, color=["#555555", *colours.values()][k], label=src)
    axes[2].axhline(0, color="k", lw=0.8)
    axes[2].set_xticks(x, groups)
    axes[2].set_ylabel("reactive minus homeostatic, mean log expression")
    st = table.set_index("source")["stress_score_mean"]
    fc = table.set_index("source")["logFC_remodelling"]
    drop = 100 * (1 - fc["harsh, high-shear"] / fc["tissue (truth)"])
    axes[2].set_title(f"Remodelling looks {drop:.0f}% weaker after harsh isolation;\nstress genes rise in every cell (mean {st['gentle, low-shear']:.2f} gentle vs {st['harsh, high-shear']:.2f} harsh)")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
