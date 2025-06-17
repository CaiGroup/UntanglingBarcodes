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

<b>readout_channels</b> - set to a python list of laser wavelengths in your experiment as annotated in the ome.tiff metadata

<b>alignment_channel</b> - set to laser wavelength of the channel to use for cross-correlation alignment (usually 405 for a dapi stain)

SeqFISH runs often probe targets sequentially (smFISH) for some hybridizations, use multiplexed combinatorial encoding for other hybridizations (seqFISH), and may use different encoding schemes for different blocks of hybridizations.

<b>initial_hyb</b> - set to the initial hybridization that encodes the seqFISH portion of your experiment.

<b>final_hyb</b> - set to the final hybridization that encodes the seqFISH portion of your experiment.

<b>z_slice</b> - set the index of the z slice to decode for (future versions will allow processing multiple z slices)

<b>positions</b> - set to "all" if you want to process all positions, or to array of non-negative integer positions in your experiment if you only want to process those positions.

## Negative Control Codewords

We add negative control codewords to our codebook that do not represent genes that we code for, so that we can use the frequency in which they are found to estimate the false postivive rate.  

<b>n_neg_cntrl_cws</b> - To use all unused codewords as negative controls, set to "all". Otherwise, set to a positive integer number less than the number of unused codewords. Setting an integer larger than or equal to the number of unused codewords is the same as setting "all".

## Preprocessing Parameters

<b>rb_radius</b> - the rolling ball radius for rolling ball background subtraction

<b>r_med_filt</b> - radius of median filter kernel applied before estimating the background with the rolling ball algorithm

<img src="https://render.githubusercontent.com/render/math?math=\left(1 - \frac{\text{hyb_num}-\text{initial_hyb}%2B1}{\text{final_hyb}-\text{initial_hyb}%2B1}\right) I %2B\left(\frac{\text{hyb_num}-\text{initial_hyb}%2B1}{\text{final_hyb} -\text{initial_hyb}%2B1}\right) F">

`white_tophat_radius`: Sometimes cells have large autofluorescent structures that drift independently of the field of view between hybridizations. You can remove thise with a white tophat of the given radius. If you do not want to use a white tophat operation, set to 'None'

## [ADCG] Parameters

<b>max_iters</b> - the maximum number of dots expected in any 64x64 pixel tile
<b>max_cd_iters</b> - the number of gradient descents to update coordinates of all dots in the model after adding each dot to the model

For the sigma_ub, sigma_lb, min_weight parameters, and final_loss_improvement use the nested entry structure described above using sigma_ub, sigma_lb, min_weight, and final_loss_improvement as super keys and readout wavelengths as subkeys. The values of these parameters must be floats (decimal) numbers. In other words, write "2.0" instead of "2".

For example: 
<pre>
sigma_ub:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 640: 2.0
</pre>


  - `sigma_ub` : for each readout_channel specified above, specify an upper bound for the psf sigma (width) fit parameter of each dot.
  - `sigma_lb` : for each readout_channel specified above, specify an lower bound for the psf sigma (width) fit parameter of each dot.
  - `min_weight` : for each readout_channel specified above, specify an lower bound for the psf sigma (width) fit parameter of each dot.
  - `min_allowed_separation` : when removing duplicates. remove dots that are within this distance of each other.


## [Syndrome Decoding] Parameters

### Cost Function Parameters

Different experiments may demand different penalty coefficients in the cost function. Set some values of the coefficients to scan. 

<b>sfs</b> - array of sigma variance penalty factors to try

<b>wfs</b> - array of log weight variance penalty factors to try

<b>lfs</b> - array of lateral variance penalty factors to try

<b>drops</b> - array of number of drops to allow: only 0 and 1 are allowed, 0 is recommended

<b>skip_thresh</b> - specifies how many barcode candidates a connected component may have before being thrown out. If "auto", set to the size of the codebook (including both on-target and off-target codewords).

## Plotting Parameters


<b>overlay_roi_width</b> - width of interative overlay to output in pixels. 200 is recommended. Large widths will produce errors.

<b>radial_density_func_rmax</b> - the maximum radius for which to calulate/plot the radial density function.

<b>radial_density_func_delta_r</b> - the delta r to use in calculating the radial density function.


[Syndrome Decoding]: https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding
[ADCG]: https://doi.org/10.1137/15M1035793