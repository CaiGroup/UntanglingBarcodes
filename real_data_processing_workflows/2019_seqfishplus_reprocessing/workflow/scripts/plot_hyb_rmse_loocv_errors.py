# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 11:12:41 2021

@author: jonat
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import numpy as np

ch = 0
wavelengths = [657, 561, 488]

for ch, ch_errors_fname in enumerate(snakemake.input):

    errors = pd.read_csv(ch_errors_fname)

    errors.loc[:,'xdiffsq'] = errors.x_diff**2
    errors.loc[:,'ydiffsq'] = errors.y_diff**2

    mean_errors = errors.groupby('hyb').mean()
    plt.plot(np.sqrt(mean_errors.xdiffsq + mean_errors.ydiffsq), label="channel %d" % wavelengths[ch])

plt.xlabel("Hybridiztion")
plt.ylabel("RMS LOOCV Alignment Error (pixels)")
plt.legend()

plt.tight_layout()
plt.savefig(snakemake.output[0])
