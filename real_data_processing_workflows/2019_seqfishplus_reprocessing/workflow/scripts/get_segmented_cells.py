# -*- coding: utf-8 -*-
"""
Created on Thu May 20 09:13:02 2021

@author: jonat
"""

import pandas as pd
import numpy as np
from skimage.draw import polygon

def get_segmented_cells(synd_dcd, pnts, rois):
    #synd_dcd['x_mean'] = [int(np.mean(eval(xs)))-1 for xs in synd_dcd['xs']]
    #synd_dcd['y_mean'] = [int(np.mean(eval(ys)))-1 for ys in synd_dcd['ys']]
    x_means = [np.mean(pnts.loc[np.array(eval(row['cpath']))-1, "x"]) for ind, row in synd_dcd.iterrows()]
    y_means = [np.mean(pnts.loc[np.array(eval(row['cpath']))-1, "y"]) for ind, row in synd_dcd.iterrows()]

    synd_dcd.loc[:, 'x_mean'] =  np.array(x_means) -1
    synd_dcd.loc[:, 'y_mean'] = np.array(y_means) -1

    synd_dcd['cell'] = -1
    
    for cell_num, cell_id in enumerate(rois.keys()):
        print("cell number: ", cell_num)
        mask = np.zeros([2048, 2048])
        _rr, _cc = polygon(rois[cell_id]['y'], rois[cell_id]['x'])
        rr = []
        cc =[]
        for i in range(len(_rr)):
            if _rr[i] < 2048 and _cc[i] < 2048:
                rr.append(_rr[i])
                cc.append(_cc[i])
        mask[rr,cc] = 1
        
        for j, barcode in synd_dcd.iterrows():
            #if mask[int(np.round(barcode['y_mean'])), int(np.round(barcode['x_mean']))]:
            x = barcode['x_mean']
            y = barcode['y_mean']
            #if barcode['y_mean'] < 2048 and barcode['x_mean'] < 2048 and mask[barcode['y_mean'], barcode['x_mean']]:
            #print("x: ", x)
            #print("y: ", y)
            if y < 2048 and x < 2048 and mask[int(y), int(x)]:
                synd_dcd.loc[j,'cell'] = cell_num
                
        sdcs = synd_dcd.loc[synd_dcd['cell'] != -1]
                
    return sdcs