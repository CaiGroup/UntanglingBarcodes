"""
Wrapper for Graph-ISS decoding (Partel & Wahlby, ICPR 2021).
Uses the author's runMaxFlowMinCost function from:
https://github.com/wahlby-lab/graph-iss/blob/master/pgm_pipeline/basecalling/pgm/max_flow_min_cost.py
"""
import numpy as np
import pandas as pd
from PIL import Image
from skimage.feature import peak_local_max
import networkx as nx
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "graph_iss_lib"))
from basecalling.pgm.max_flow_min_cost import runMaxFlowMinCost

sys.path.insert(0, os.path.dirname(__file__))
from convert_codebook import load_codebook

cb_path = snakemake.input[0]
code = snakemake.params["code"]
zeros_probed = code[:7] == "seqFISH"
cb, matrix, n_rounds, q, gene_col = load_codebook(cb_path)
n_channels = q if zeros_probed else max(q - 1, 1)

# seqFISH-family codebooks have no blank symbol: every round always carries a real
# dot, and symbol 0 is rendered/detected as pseudocolor q rather than being skipped
# (see gen_sim_data.jl and convert_codebook.py). Match against the same remapped
# symbols the detected spots use, so a detected ch=q correctly matches a codebook
# entry of 0.
match_matrix = np.where(matrix == 0, q, matrix) if zeros_probed else matrix

image_files = snakemake.input[1:]
first_img = np.array(Image.open(image_files[0]), dtype=np.float64)
H, W = first_img.shape

images = {}
for img_path in image_files:
    fname = os.path.basename(img_path)
    parts = fname.split("_")
    r = int(parts[3])
    pc = int(parts[5])
    img = np.array(Image.open(img_path), dtype=np.float64)
    img_max = img.max()
    if img_max > 0:
        img = img / img_max
    images[(r, pc)] = img

min_distance = snakemake.params["spot_min_distance"]
threshold_rel = snakemake.params["spot_threshold_rel"]
spatial_radius = snakemake.params["spatial_radius"]

spots = []
for (r, pc), img in images.items():
    if img.max() == 0:
        continue
    thresh = img.max() * threshold_rel
    coords = peak_local_max(img, min_distance=min_distance, threshold_abs=thresh)
    for coord in coords:
        y_pos, x_pos = coord
        spots.append({
            "hyb": r,
            "ch": pc,
            "x": float(x_pos),
            "y": float(y_pos),
            "z": 0.0,
            "p0": 0.1,
            "p1": 0.9,
            "Imax_gf": float(img[y_pos, x_pos])
        })

if len(spots) == 0:
    pd.DataFrame(columns=["gene", "x", "y"]).to_csv(snakemake.output[0], index=False)
else:
    data = pd.DataFrame(spots)

    tagList = pd.DataFrame(match_matrix.astype(int))

    d_th = spatial_radius
    dth_max = spatial_radius
    k1 = snakemake.params["k1"]

    res = runMaxFlowMinCost(data, d_th, k1, tagList, n_threads=1, dth_max=dth_max, prior="prior")

    decoded_genes = []
    decoded_xs = []
    decoded_ys = []

    if res is not None:
        num_hybs = n_rounds
        for component in res:
            G = component['G']
            Dvar = component['Dvar']
            for cc in nx.connected_components(G):
                cc_arr = np.array(list(cc))
                cc_arr = cc_arr[cc_arr <= Dvar.X_idx.max()]
                dvar_cc = Dvar[Dvar.X_idx.isin(cc_arr)]
                if len(dvar_cc) == num_hybs:
                    barcode = [0] * num_hybs
                    for _, row in dvar_cc.iterrows():
                        barcode[int(row.hyb) - 1] = int(row.ch)
                    barcode_tuple = tuple(barcode)
                    for j in range(len(match_matrix)):
                        if tuple(int(match_matrix[j, r]) for r in range(n_rounds)) == barcode_tuple:
                            decoded_genes.append(int(cb.iloc[j][gene_col]))
                            decoded_xs.append(float(dvar_cc.iloc[0].x))
                            decoded_ys.append(float(dvar_cc.iloc[0].y))
                            break

    results = pd.DataFrame({
        "gene": decoded_genes,
        "x": decoded_xs,
        "y": decoded_ys
    })
    results.to_csv(snakemake.output[0], index=False)
