import pandas as pd
import numpy as np
import os
import re

benchmark_dirs = {
    "Untangle: spots first": [
        "results/benchmarks/spotsfirst_candidates",
        "results/benchmarks/spotsfirst_ilp"
    ],
    "Untangle: CS l0": [
        "results/benchmarks/cs_spot_detection",
        "results/benchmarks/cs_cand_cpaths",
        "results/benchmarks/cs_cand_cpaths_nc",
        "results/benchmarks/cs_l0learn"
    ],
    "bardensr": ["results/benchmarks/bardensr"],
    "istdeco": ["results/benchmarks/istdeco"],
    "graph_iss": ["results/benchmarks/graph_iss"],
    "polaris": ["results/benchmarks/polaris"],
}

rows = []
for method, dirs in benchmark_dirs.items():
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if not f.endswith(".tsv"):
                continue
            df = pd.read_csv(os.path.join(d, f), sep="\t")
            m = re.match(
                r"(.+?)_rep_(\d+)_nbarcodes(\d+)_nnonspec(\d+)_pdrop([\d.]+)_rstdv([\d.]+)",
                f.replace(".tsv", "")
            )
            if m:
                code = m.group(1)
                rep = int(m.group(2))
                nbarcodes = int(m.group(3))
                rstdv = float(m.group(6))
            else:
                code, rep, nbarcodes, rstdv = f, 0, 0, 0.0
            rows.append({
                "method": method,
                "step": os.path.basename(d),
                "code": code,
                "rep": rep,
                "nbarcodes": nbarcodes,
                "rstdv": rstdv,
                "seconds": df["s"].iloc[0],
                "max_rss_mb": df["max_rss"].iloc[0] if "max_rss" in df.columns else np.nan,
                "cpu_time": df["cpu_time"].iloc[0] if "cpu_time" in df.columns else np.nan,
            })

result = pd.DataFrame(rows)

total_time = result.groupby(
    ["method", "code", "rep", "nbarcodes", "rstdv"]
).agg(
    total_seconds=("seconds", "sum"),
    max_rss_mb=("max_rss_mb", "max"),
).reset_index()

total_time.to_csv(snakemake.output[0], index=False)
