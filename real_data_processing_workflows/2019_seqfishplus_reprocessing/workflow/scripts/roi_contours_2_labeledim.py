from read_roi import read_roi_zip
import numpy as np
from skimage.draw import polygon
import tifffile


rois_filename = snakemake.input[0]#"resources/rep1_Roi/RoiSet_Pos0.zip"
rois = read_roi_zip(rois_filename)

labeled_im = np.zeros([2048, 2048], dtype=np.uint8)

for cell_num, cell_id in enumerate(rois.keys()):
    print("cell number: ", cell_num)
    _rr, _cc = polygon(rois[cell_id]['y'], rois[cell_id]['x'])
    rr = []
    cc =[]
    for i in range(len(_rr)):
        if _rr[i] < 2048 and _cc[i] < 2048:
            rr.append(_rr[i])
            cc.append(_cc[i])
    labeled_im[rr,cc] = cell_num + 1

tifffile.imsave(snakemake.output[0], np.uint8(labeled_im))
