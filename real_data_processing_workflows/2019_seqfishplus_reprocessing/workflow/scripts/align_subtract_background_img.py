
import numpy as np
import pandas as pd

"""
align_subtract_background_img(hyb_im, bg, offsets, hyb)

Aligns and subtracts an initial background and fiducial marker image from an image of a readout hybridization.
"""
def align_subtract_background_img(hyb_im, bg, offsets, hyb):
    x_t = offsets.loc[hyb, "x"]
    y_t = offsets.loc[hyb, "y"]

    x_p = x_t % 1
    y_p = y_t % 1

    x_tf = int(np.floor(x_t))
    x_tc = int(np.ceil(x_t))
    y_tf = int(np.floor(y_t))
    y_tc = int(np.ceil(y_t))

    xsize, ysize = np.shape(bg)
    bgt = np.zeros((xsize, ysize))

    for i in range(xsize):
        for j in range(ysize):
            if j - x_tc < xsize and i - y_tc < ysize and j - x_tc > 0 and i - y_tc > 0 and j - x_tf < xsize and i - y_tf < ysize and j - x_tf > 0 and i - y_tf > 0:
                bgt[i, j] += bg[i - y_tc, j - x_tc] * x_p * y_p
                bgt[i, j] += bg[i - y_tf, j - x_tc] * x_p * (1 - y_p)
                bgt[i, j] += bg[i - y_tc, j - x_tf] * (1 - x_p) * y_p
                bgt[i, j] += bg[i - y_tf, j - x_tf] * (1 - x_p) * (1 - y_p)

    im_bgt_sub = hyb_im - bgt
    im_bgt_sub[im_bgt_sub < 0] = 0

    return im_bgt_sub

def translate_img(img, offsets, hyb):
    x_t = offsets.loc[hyb, "x"]
    y_t = offsets.loc[hyb, "y"]

    x_p = x_t % 1
    y_p = y_t % 1

    x_tf = int(np.floor(x_t))
    x_tc = int(np.ceil(x_t))
    y_tf = int(np.floor(y_t))
    y_tc = int(np.ceil(y_t))

    xsize, ysize = np.shape(img)
    translated = np.zeros((xsize, ysize))

    for i in range(xsize):
        for j in range(ysize):
            if j - x_tc < xsize and i - y_tc < ysize and j - x_tc > 0 and i - y_tc > 0 and j - x_tf < xsize and i - y_tf < ysize and j - x_tf > 0 and i - y_tf > 0:
                translated[i, j] += img[i - y_tc, j - x_tc] * x_p * y_p
                translated[i, j] += img[i - y_tf, j - x_tc] * x_p * (1 - y_p)
                translated[i, j] += img[i - y_tc, j - x_tf] * (1 - x_p) * y_p
                translated[i, j] += img[i - y_tf, j - x_tf] * (1 - x_p) * (1 - y_p)

    return translated