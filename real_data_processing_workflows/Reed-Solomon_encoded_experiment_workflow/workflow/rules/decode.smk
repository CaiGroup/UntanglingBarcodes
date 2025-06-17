import skimage
import numpy as np

for pos in config["positions"]:
    #rule split_cells
    lbld_img = skimage.io.imread("resources/labeled_images/pos{pos}_polyT_dapi_cp_masks.png".format(pos=pos))
    np.unique(lbld_img)
    cells = [cellID for cellID in np.unique(lbld_img) if cellID != 0]
    rule:
        input: "results/alignment/ref_cross_channel_registered_dots_pos{pos}.csv".format(pos=pos)
        output: expand("results/cell_dots/cell_dots_pos"+str(pos)+"_cell{{cell}}.csv", cell=cells)
        resources: mem_mb = 5000
        script: "../scripts/split_cells.py"

    #agregate decoded_cell
    rule:
        input: expand("results/decoded/barcodes_unfiltered_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells), 
	wildcard_constraints: cell="\d+", pos="\d+", dr="\d"
        output: "results/decoded/barcodes_unfiltered_pos_{pos}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos)
        resources: mem_mb = 5000
        script: "../scripts/aggregate_sum_stats.py"

    rule:
        input: expand("results/sum_stats/sum_stats_unfiltered_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells), 
        output: "results/sum_stats/sum_stats_unfiltered_pos_{pos}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos)
	wildcard_constraints: pos="\d+"
        resources: mem_mb = 5000
        script: "../scripts/aggregate_sum_stats.py"


rule find_codepaths:
    input: 
           "results/cell_dots/cell_dots_pos{pos}_cell{cell}.csv",
           "resources/codebooks/full_RS_q11_k7_cb.csv",
           "resources/codebooks/RS_q11_k7_H.csv",
    #params: lf=min(lfs), wf = min(wfs), sf = min(sfs), dr = max(drops)
    params: lf = min(config['lfs']),
            zf = min(config['zfs']),
            wf = min(config['wfs']),
            sf = min(config['sfs']),
            dr = max(config['drops']),
            rxy = config["rxy_ro"], 
            rz = config["rz_ro"]
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    output: "results/codepaths/cps_pos_{pos}_cell{cell}.csv"
    group: "find codepaths"
    resources: mem_mb=20000
    script: "../scripts/find_save_codepath_candidates.jl"


rule choose_optimal_codepaths_w_neg_ctrl:
    input:
        "results/codepaths/cps_pos_{pos}_cell{cell}.csv",
        "results/cell_dots/cell_dots_pos{pos}_cell{cell}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv",
        "resources/codebooks/RS_q11_k7_H.csv",
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config["rxy_ro"], rz = config["rz_ro"]
    group: "decode"
    output:
        "results/decoded/barcodes_unfiltered_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_w_neg_ctrl_unfiltered_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/decoded/points_pos_{pos}_cell{cell}_lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/dense_discarded_cpaths/dense_cpaths__pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
        #"results/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        #temp("results/sum_stats/sum_stats_unfiltered_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"),
        #"results/decoded/points_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        #"results/dense_discarded_cpaths/dense_cpaths_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
        #"results/rep{rep}/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        #temp("results/rep{rep}/sum_stats/sum_stats_unfiltered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"),
        #temp("results/rep{rep}/decoded/points_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv")
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", cell="\d+"
    resources: mem_mb = 40000
    script: "../scripts/choose_cpaths_from_saved_candidates.jl"

rule choose_optimal_codepaths_no_neg_ctrl:
    input:
        "results/codepaths/cps_pos_{pos}_cell{cell}.csv",
        "results/cell_dots/cell_dots_pos{pos}_cell{cell}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv",
        "resources/codebooks/RS_q11_k7_H.csv",
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config["rxy_ro"], rz = config["rz_ro"]
    output:
        "results/decoded/barcodes_unfiltered_no_neg_ctrl_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_unfiltered_no_neg_ctrl_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/decoded/points_no_neg_ctrl_pos_{pos}_cell{cell}_lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/dense_discarded_cpaths/dense_cpaths_no_neg_ctrl_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", cell="\d+"
    resources: mem_mb = 20000
    script: "../scripts/choose_cpaths_from_saved_candidates_no_neg_ctrl.jl"

rule reconcile_sum_stats:
    input:
        "results/sum_stats/sum_stats_w_neg_ctrl_unfiltered_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_unfiltered_no_neg_ctrl_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv"
    output: "results/sum_stats/sum_stats_unfiltered_pos_{pos}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    resources: mem_mb=1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/estimate_fdr.py"


rule aggregate_sum_stats:
    #input: expand("results/sum_stats/sum_stats_unfiltered_pos_{pos}_mw_{mw}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv", lf=config["lfs"], zf=config["zfs"], wf=config["wfs"], sf=config["sfs"], dr=config["drops"], pos=config["positions"], ch=config["channels"], mw=config["min_weights"])
    input: expand("results/sum_stats/sum_stats_unfiltered_pos_{pos}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv", lf=config["lfs"], zf=config["zfs"], wf=config["wfs"], sf=config["sfs"], dr=config["drops"], pos=config["positions"], ch=config["channels"])
    #expand("results/sum_stats/sum_stats_pos_{pos}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv", lf = lfs, wf = wfs, sf = sfs, dr = drops, pos=positions, zf=zfs)
    output: "results/decode_summary_stats_pos_{pos}.csv"
    script: "../scripts/aggregate_sum_stats.py"


rule aggregate_decoded_position:
    input:
        "resources/pos.pos",
        expand('results/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv', ch=config['channels'], pos=config['positions'])
    output:
        "results/decoded_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    script: "../scripts/aggregate_positions.py"




"""
rule apply_xcorr_offsets:
    input: 'alignments/xcorr_alignments.csv', "fit_dots/all_hyb_dots.csv"
    output: "alignments/aligned_dots.csv"
    script: "scripts/apply_xcorr_alignments.py"
"""

"""
rule aggregate_decoded_position:
    input:
        "resources/pos.pos",
        expand('results/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv', ch=config['channels'], pos=config['positions'])
    output:
        "results/decoded_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    script: "../scripts/aggregate_positions.py"
"""

rule filter_pos_ch_barcodes:
    input:
        "results/rep{rep}/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv"
    params: filt_size=config["min_filter_size"], filt_prop=config["filter_prop_neg_control"]
    group: "decode"
    output:
        "results/rep{rep}/decoded/barcodes_filtered_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/sum_stats/sum_stats_filtered_ch_{ch}_pos_{pos}_mw_{mw}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb=1000
    script: "../scripts/filter_barcodes.py"
