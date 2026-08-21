import numpy as np
import pandas as pd
from PIL import Image
import sys
import os
from skimage.feature import peak_local_max

sys.path.insert(0, os.path.dirname(__file__))
from convert_codebook import to_istdeco, load_codebook

istdeco_path = snakemake.params["istdeco_path"]
if not os.path.isabs(istdeco_path):
    istdeco_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(istdeco_path))
sys.path.insert(0, istdeco_path)
from istdeco import ISTDeco

cb_path = snakemake.input[0]
cb, matrix, n_rounds, q, gene_col = load_codebook(cb_path)
n_channels = max(q - 1, 1)

image_files = snakemake.input[1:]
first_img = np.array(Image.open(image_files[0]), dtype=np.float32)
H, W = first_img.shape

Y = np.zeros((n_rounds, n_channels, H, W), dtype=np.float32)

for img_path in image_files:
    fname = os.path.basename(img_path)
    parts = fname.split("_")
    r = int(parts[3])
    pc = int(parts[5])
    img = np.array(Image.open(img_path), dtype=np.float32)
    if q == 2:
        Y[r - 1, 0, :, :] = img
    else:
        Y[r - 1, pc - 1, :, :] = img

y_max = Y.max()
if y_max > 0:
    Y = Y / y_max

D = to_istdeco(cb_path)
D = D.astype(np.float32)

sigma = snakemake.params["sigma"]
b = snakemake.params["b"]
niter = snakemake.params["niter"]
suppress_radius = snakemake.params["suppress_radius"]
device = snakemake.params["device"]
combos = snakemake.params["combos"]
outputs = snakemake.output

import torch
model = ISTDeco(Y, D, (sigma, sigma), b=b)
model = model.to(device)
X, Q, loss = model.run(niter=niter, suppress_radius=suppress_radius)

nonzero_Y = Y[Y > 0]

for idx, (threshold_Q, threshold_percentile) in enumerate(combos):
    tau = np.percentile(nonzero_Y, threshold_percentile) if len(nonzero_Y) > 0 else 0.01

    genes = []
    xs = []
    ys = []
    intensities = []
    qualities = []

    for j in range(X.shape[0]):
        mask = (X[j] > tau) & (Q[j] > threshold_Q)
        if not np.any(mask):
            continue
        coords = peak_local_max(X[j] * mask.astype(float), min_distance=1, threshold_abs=tau)
        for coord in coords:
            y_pos, x_pos = coord
            genes.append(int(cb.iloc[j][gene_col]))
            xs.append(float(x_pos))
            ys.append(float(y_pos))
            intensities.append(float(X[j, y_pos, x_pos]))
            qualities.append(float(Q[j, y_pos, x_pos]))

    results = pd.DataFrame({
        "gene": genes,
        "x": xs,
        "y": ys,
        "intensity": intensities,
        "quality": qualities
    })

    results.to_csv(outputs[idx], index=False)
