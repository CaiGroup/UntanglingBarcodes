# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 08:15:58 2021

@author: jonat
"""

from scipy.linalg import circulant
import numpy as np
from skimage.draw import polygon2mask
from skimage.io import imread

def auto_roi(aligned_hybs, fov_size, rois=None):

    stack_sum = np.sum(aligned_hybs, axis=0)

    xdim, ydim = np.shape(stack_sum)

    if rois:
        mask = np.zeros((ydim, xdim))

        for cell_num, cell_id in enumerate(rois.keys()):
            cell_dict = rois[cell_id]
            contour = tuple(zip(cell_dict['y'], cell_dict['x']))
            mask += polygon2mask((ydim, xdim), contour)

        mask[mask>1] = 1

        masked_stack_sum = mask*stack_sum
    else:
        masked_stack_sum = stack_sum

    # draw roi mask

    # get circulat matrices
    c1 = np.zeros(2048)
    c1[:fov_size] = 1

    to_trim = fov_size-1

    m = circulant(c1)[to_trim:]
    mt = np.transpose(m)

    ul_corner_overlap = np.dot(np.dot(m,masked_stack_sum),mt)

    max_overlaps = np.argwhere(ul_corner_overlap == np.max(ul_corner_overlap))

    return max_overlaps[0]

    #for cell_num, cell_id in enumerate(rois.keys()):
