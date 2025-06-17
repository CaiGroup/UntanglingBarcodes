from util import pil_imread
import xarray
import numpy as np
import tifffile
import skimage
import skimage.io

tifstack_name = snakemake.input[0]

if "z" in snakemake.wildcards.keys():
    z = int(snakemake.wildcards["z"])
else:
    z = snakemake.params["z"]

if "ch" in snakemake.wildcards.keys():
    save_ch = int(snakemake.wildcards["ch"])
else:
    save_ch = snakemake.params["ch"]

outfname = snakemake.output[0]

tifstack, metadata = pil_imread(tifstack_name,metadata=True)

channels = np.unique([int(immddict['Channel']) for immddict in metadata])
channels=-np.sort(-channels)
nchannels, nzs, nys, nxs = np.shape(tifstack)
ds = xarray.Dataset({'fluorescence': (("channel", "z", "y", "x"), tifstack)},
                         {"channel":channels,
                          "z": range(nzs),
                          "y": range(nys),
                          "x": range(nxs)
                         })

#max_z_projection = ds.max(dim="z")
#ch_im = max_z_projection.sel(channel=save_ch)
#max_z_projection = ds.max(dim="z")
#ch_im = max_z_projection.sel(channel=save_ch)
ch_im = ds.sel(channel=save_ch, z=z)

#tifffile.imsave(outfname, ch_im.fluorescence)
skimage.io.imsave(outfname, ch_im.fluorescence)
