import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

data = pd.read_csv(snakemake.input[0])
data = data.dropna(subset=["method"])

methods = sorted(data["method"].unique())
densities = sorted(data["nbarcodes"].unique())
codes = sorted(data["code"].unique())

colors = plt.cm.tab10(np.linspace(0, 1, max(len(methods), 10)))
method_colors = {m: colors[i] for i, m in enumerate(methods)}

markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":"]
density_markers = {d: markers[i % len(markers)] for i, d in enumerate(densities)}
density_linestyles = {d: linestyles[i % len(linestyles)] for i, d in enumerate(densities)}

sweep_col = "rstdv"
sweep_label = "Localization jitter (std dev)"

n_codes = len(codes)
fig, axes = plt.subplots(1, n_codes, figsize=(7 * n_codes, 6), squeeze=False)

for ci, code in enumerate(codes):
    ax = axes[0, ci]
    cdata = data[data["code"] == code]
    for method in methods:
        for density in densities:
            mdata = cdata[(cdata["method"] == method) & (cdata["nbarcodes"] == density)]
            if len(mdata) == 0:
                continue
            grouped = mdata.groupby(sweep_col).agg(
                time_mean=("total_seconds", "mean"),
                time_sem=("total_seconds", "sem"),
            ).reset_index()
            label = f"{method} (n={int(density)})"
            ax.errorbar(grouped[sweep_col], grouped["time_mean"],
                         yerr=grouped["time_sem"],
                         label=label, color=method_colors[method],
                         marker=density_markers[density],
                         linestyle=density_linestyles[density],
                         capsize=3, markersize=5)
    ax.set_xlabel(sweep_label)
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title(code)
    ax.legend(fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150)
plt.savefig(snakemake.output[1])
