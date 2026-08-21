import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

data = pd.read_csv(snakemake.input[0])
data = data.dropna(subset=["method"])
code = snakemake.wildcards["code"]

methods = sorted(data["method"].unique())
densities = sorted(data["nbarcodes"].unique())

colors = plt.cm.tab10(np.linspace(0, 1, max(len(methods), 10)))
method_colors = {m: colors[i] for i, m in enumerate(methods)}

markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":"]
density_markers = {d: markers[i % len(markers)] for i, d in enumerate(densities)}
density_linestyles = {d: linestyles[i % len(linestyles)] for i, d in enumerate(densities)}

sweep_col = "rstdv"
sweep_label = "Localization jitter (std dev)"

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for method in methods:
    for density in densities:
        mdata = data[(data["method"] == method) & (data["nbarcodes"] == density)]
        if len(mdata) == 0:
            continue
        grouped = mdata.groupby(sweep_col).agg(
            sens_mean=("sensitivity", "mean"),
            sens_sem=("sensitivity", "sem"),
            fdr_mean=("fdr", "mean"),
            fdr_sem=("fdr", "sem")
        ).reset_index()

        x = grouped[sweep_col].values
        label = f"{method} (n={int(density)})"
        axes[0].errorbar(x, grouped["sens_mean"], yerr=grouped["sens_sem"],
                         label=label, color=method_colors[method],
                         marker=density_markers[density],
                         linestyle=density_linestyles[density],
                         capsize=3, markersize=5)
        axes[1].errorbar(x, grouped["fdr_mean"], yerr=grouped["fdr_sem"],
                         label=label, color=method_colors[method],
                         marker=density_markers[density],
                         linestyle=density_linestyles[density],
                         capsize=3, markersize=5)

axes[0].set_xlabel(sweep_label)
axes[0].set_ylabel("Sensitivity")
#axes[0].set_ylim(-0.05, 1.05)
axes[0].legend(fontsize=7, ncol=2, loc="lower left")
axes[0].set_title(f"Sensitivity — {code}")

axes[1].set_xlabel(sweep_label)
axes[1].set_ylabel("FDR")
#axes[1].set_ylim(-0.05, 1.05)
axes[1].legend(fontsize=7, ncol=2, loc="upper left")
axes[1].set_title(f"FDR — {code}")

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150)
plt.savefig(snakemake.output[1])
