# Output Files

When a run finishes, the workflow will save two .csv files in the results directory: <i>decoded_summary_stats.csv </i> and <i> decoded_barcode.csv </i>.

<i>decoded_summary_stats.csv </i> gives a summary of how well each position-channel image decoded with each set of syndrome decoding parameters.

There tends to be a trade off between how many on target barcodes are found and the estimated false positive rate in setting the syndrome decoding parameters. It is up to you to decide how many false positives you are willing to tolerate to get more true positive barcodes.

To you help decide which syndrome decoding parameters to use, there are plots of the the on-target/off-target counts for each syndrome decoding parameter set used for each position-channel image in the results/plots/on_off_plots folder.

Also to help gauge the performance the workflow, there will be interactive plots showing overlays of dots found in ADCG and how they were or were not decoded in the results/interactive_plots folder. Mouse over on-target dots to see what gene they were decoded as and what their codeword is in hovertext.

<i> results/mpaths_ch_*_pos_*_lf*_wf*_sf*_df*.csv </i> files list each decoded barcode in the experiment with many columns:

  * <b>gene</b> gives the name of the gene for which the barcode represents (negative control means that barcode is for codeword that was unassigned to a gene)
  * <b>gene_number</b> is the row which the genes (or negative control's) codeword was listed in the codebook
  * <b>cpath</b> is a list of the indices of the dots for the cell that they are listed in the aligned_dots_ch_*_pos_*.csv
  * <b>cost</b> is the cost that the barcode was assigned by the cost funtion during decoding
  * <b>cellid</b> is the id of the cell in that particuar position that the barcode was found
  * <b>x</b> is the average x coordinate of the aligned dots in the barcode
  * <b>y</b> is the average y coordinate of the aligned dots in the barcode
  * <b>cc</b> is the connected component (of candidate barcodes) in which the barcode was chosen during decoding
  * <b>cc_size</b> is the number of candidate barcodes that were in the connected component of candidate barcodes in which the barcode was chosen during decoding
  * <b>pos</b> is the position number in which the barcode was found
  * <b>ch</b> is channel in which the barcode was found
  * <b>lvf</b> is the lateral variance penalty factor used in the syndrome decoding cost function
  * <b>wf</b> is log weight variance penalty factor used in the syndrome decoding cost function
  * <b>svf</b> is the sigma variance penalty factor used in the syndrome decoding cost function
  * <b>dr</b> is the number of drops allowed in a decoded barcode
  
 You may also wish to remove on-target barcodes that were found in connected components from which many negative control barcodes were also found.