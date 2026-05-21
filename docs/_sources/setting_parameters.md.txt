# Setting Parameters

Set the parameters for the processing run in the [config/config.yaml] file.

.yaml files can be interpreted as specifying dictionaries of the form

<pre> key: value </pre>

entries can be nested, for example

<pre>
super_key: <n>
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; sub_key_1: value
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; sub_key_2: value
</pre>
## File organization parameters
* <b>readout_channels</b> - set to a python list of laser wavelengths in your experiment as annotated in the ome.tiff metadata
* <b>alignment_channel</b> - set to laser wavelength of the channel to use for cross-correlation alignment (usually 405 for a dapi stain)

SeqFISH runs often probe targets sequentially (smFISH) for some hybridizations, use multiplexed combinatorial encoding for other hybridizations (seqFISH), and may use different encoding schemes for different blocks of hybridizations.
* <b>initial_hyb</b> - set to the initial hybridization that encodes the seqFISH portion of your experiment.
* <b>final_hyb</b> - set to the final hybridization that encodes the seqFISH portion of your experiment.
* <b>z_slice</b> - set the index of the z slice to decode for (future versions will allow processing multiple z slices)
* <b>positions</b> - set to "all" if you want to process all positions, or to array of non-negative integer positions in your experiment if you only want to process those positions.

## Negative Control Codewords

We add negative control codewords to our codebook that do not represent genes that we code for, so that we can use the frequency in which they are found to estimate the false postivive rate.  

<b>n_neg_cntrl_cws</b> - To use all unused codewords as negative controls, set to "all". Otherwise, set to a positive integer number less than the number of unused codewords. Setting an integer larger than or equal to the number of unused codewords is the same as setting "all".

## Preprocessing Parameters

<b>rb_radius</b> - the rolling ball radius for rolling ball background subtraction

<b>r_med_filt</b> - radius of median filter kernel applied before estimating the background with the rolling ball algorithm

### Reed-Solomon Workflow only
* <b>bg_sub_multiplier</b> - the aligned background is multiplied by this number when subtracted from each readout image. It is set to 1 in final version of the workflow

Beads used for fiducial marker alignement overlap with cells in the Reed-Solomon encoded dataset and require an extra step to mask them out in order to entirely remove them from the images before fitting. The following parameters are used to find the bead masks:
* <b>bead_thresh</b> The threshold above which pixels are assumed to be bead in the bead only images. The initial bead mask includes pixels above this threshold.
* <b>bead_dilation_radius</b> - This sets a dialation radius for the bead mask increase the masking area and ensure that the beads are completely masked out.

## ADCG Parameters

Our implementation of ADCG divides images into small overlapping square tiles of width `tile_main_width + tile_overlap` to improve computational efficiency
* <b>tile_main_width</b> - The spacing between starts of tiles
* <b>tile_overlap</b> - the length of the overlap regions between adjacent tiles. It should be at least long enough to enclose an entire point spread function.
* <b>max_iters</b> - the maximum number of dots expected in any 64x64 pixel tile
* <b>max_cd_iters</b> - the number of gradient descents to update coordinates of all dots in the model after adding each dot to the model

For the sigma_ub, sigma_lb, min_weight parameters, and final_loss_improvement use the nested entry structure described above using sigma_ub, sigma_lb, min_weight, and final_loss_improvement as super keys and readout wavelengths as subkeys. The values of these parameters must be floats (decimal) numbers. In other words, write "2.0" instead of "2".

For example: 
<pre>
sigma_ub:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 640: 2.0
</pre>


  - <b>sigma_ub</b> - for each readout_channel specified above, specify an upper bound for the psf sigma (width) fit parameter of each dot.
  - <b>sigma_lb</b> - for each readout_channel specified above, specify an lower bound for the psf sigma (width) fit parameter of each dot.
  - <b>min_weight</b> - for each readout_channel specified above, specify an lower bound for the psf sigma (width) fit parameter of each dot.
  - <b>min_allowed_separation</b> - when removing duplicates. remove dots that are within this distance of each other.

## Fiducial marker matching parameters.
* <b>min_fm_hyb_matches</b>
* <b>outlier_sd_thresh</b>
* <b>max_lat_offset</b>
* <b>set_xy_search_error</b>

## Syndrome Decoding Parameters

* <b>r_xy</b> - search radius for fiducial markers
* <b>r_z</b> - z search radius for fiducial markers (not used)
* <b>rxy_ro</b>  - search radius for decoding readout dots
* <b>rz_ro</b> - z search radius for decoding readout dots (not used)

### Cost Function Parameters

Different experiments may demand different penalty coefficients in the cost function. Set some values of the coefficients to scan. 
* <b>sfs</b> - array of sigma variance penalty factors to try
* <b>wfs</b> - array of log weight variance penalty factors to try
* <b>lfs</b> - array of lateral variance penalty factors to try
* <b>drops</b> - array of number of drops to allow: only 0 and 1 are allowed, 0 is recommended
* <b>skip_thresh</b> - specifies how many barcode candidates a connected component may have before being thrown out. If "auto", set to the size of the codebook (including both on-target and off-target codewords).

### Parameters for deciding whether or not to discard a decoding conflict network
* <b>skip_thresh</b> - Discard conflict networks containing more than skip_thresh barcodes and a higher density per pixel than skip_density_thresh
* <b>skip_density_thresh</b> - Discard conflict networks containing more than skip_thresh barcodes and a higher density per pixel than skip_density_thresh
* <b>min_filter_size</b> - A conservative approach to decoding is to use decoded genes decoded at the same time as negative control codewords. When doing this, you may choose 
* <b>filter_prop_neg_control</b> 0.3

## Plotting Parameters
* <b>overlay_roi_width</b> - width of interative overlay to output in pixels. 200 is recommended. Large widths will produce errors.
* <b>radial_density_func_rmax</b> - the maximum radius for which to calulate/plot the radial density function.
* <b>radial_density_func_delta_r</b> - the delta r to use in calculating the radial density function.


[Syndrome Decoding]: https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding
[ADCG]: https://doi.org/10.1137/15M1035793