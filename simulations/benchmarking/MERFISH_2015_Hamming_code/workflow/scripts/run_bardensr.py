import numpy as np
import pandas as pd
from PIL import Image
import bardensr.spot_calling as bardensr
import sys
import os
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from convert_codebook import to_bardensr, load_codebook

cb_path = snakemake.input[0]
cb, matrix, n_rounds, q, gene_col = load_codebook(cb_path)
n_channels = max(q - 1, 1)
n_images = n_rounds * n_channels if q > 2 else n_rounds

image_files = snakemake.input[1:]

first_img = np.array(Image.open(image_files[0]), dtype=np.float64)
H, W = first_img.shape

imagestack = np.zeros((n_images, 1, H, W), dtype=np.float64)

for img_path in image_files:
    fname = os.path.basename(img_path)
    parts = fname.split("_")
    r = int(parts[3])
    pc = int(parts[5])
    if q == 2:
        frame_idx = r - 1
    else:
        frame_idx = (r - 1) * n_channels + (pc - 1)
    img = np.array(Image.open(img_path), dtype=np.float64)
    imagestack[frame_idx, 0, :, :] = img

img_max = imagestack.max()
if img_max > 0:
    imagestack = imagestack / img_max

codebook = to_bardensr(cb_path)

sigma = snakemake.params["sigma"]
iterations = snakemake.params["iterations"]
poolsize = snakemake.params["poolsize"]
combos = snakemake.params["combos"]
outputs = snakemake.output

# Group (peak_thresh, output path) by l1 so estimate_density_iterative -- the expensive step --
# runs once per distinct l1 penalty instead of once per (l1, peak) combo.
peaks_by_l1 = defaultdict(list)
for idx, (l1_penalty, peak_thresh) in enumerate(combos):
    peaks_by_l1[l1_penalty].append((peak_thresh, outputs[idx]))

for l1_penalty, peak_outputs in peaks_by_l1.items():
    density, gains = bardensr.estimate_density_iterative(
        imagestack, codebook,
        l1_penalty=l1_penalty,
        psf_radius=(0, int(round(sigma)), int(round(sigma))),
        iterations=iterations,
        estimate_codebook_gain=True
    )

    for peak_thresh, out_path in peak_outputs:
        if peak_thresh == "auto":
            nonzero_vals = density[density > 0]
            if len(nonzero_vals) > 0:
                thresh = np.percentile(nonzero_vals, 90)
            else:
                thresh = 0.01
        else:
            thresh = float(peak_thresh)

        spots_df = bardensr.find_peaks(
            density, thresh,
            poolsize=(1, poolsize, poolsize)
        )

        if len(spots_df) > 0:
            magnitudes = [float(density[int(row.m0), int(row.m1), int(row.m2), int(row.j)]) for _, row in spots_df.iterrows()]
            results = pd.DataFrame({
                "gene": spots_df["j"].values + 1,
                "x": spots_df["m2"].values.astype(float),
                "y": spots_df["m1"].values.astype(float),
                "magnitude": magnitudes
            })
        else:
            results = pd.DataFrame(columns=["gene", "x", "y", "magnitude"])

        results.to_csv(out_path, index=False)
