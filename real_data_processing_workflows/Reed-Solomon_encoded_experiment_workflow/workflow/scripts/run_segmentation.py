import tifffile
from cellpose import models
import util
import numpy as np

#img = skimage.io.imread("MMStacks/DMSO2/HybCycle_40/MMStack_Pos%d.ome.tif" % pos)

print("loading stack")
print(snakemake.input[0])
hyb_stack, md = util.pil_imread(snakemake.input[0], metadata=True)

print(md)
ind_2_ch_dict = {}
for im_md in md: ind_2_ch_dict[im_md["ChannelIndex"]] = im_md["ChannelIndex"] 

print(ind_2_ch_dict)
print(np.shape(hyb_stack))
poly_t_stack = hyb_stack[2,2,:,:]
    
model = models.Cellpose(model_type='cyto2')

print("running model")
results = model.eval(poly_t_stack, diameter=300, channels=[[0,0]], z_axis=0)
labeled_stack = results[0]

tifffile.imwrite(snakemake.output[0], labeled_stack)