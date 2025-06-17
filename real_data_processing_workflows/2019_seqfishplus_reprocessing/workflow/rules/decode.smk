import skimage
import numpy as np



for rep in config['positions']:
    for pos in config["positions"][rep]:
        #rule split_cells
        lbld_img = skimage.io.imread("resources/replicates/"+rep+"/Labeled_Images/MMStack_Pos{pos}.png".format(pos=pos))
        np.unique(lbld_img)
        cells = [cellID for cellID in np.unique(lbld_img) if cellID != 0]
        rule:
            input: "results/{rep}/alignment/aligned_dots_ch_{{ch}}_pos_{pos}.csv".format(pos=pos, rep=rep)
            output: expand("results/"+rep+"/cell_dots/cell_dots_ch_{{ch}}_pos_"+str(pos)+"_cell{cell}.csv", cell=cells)
            resources: mem_mb = 5000
            script: "../scripts/split_cells.py"

        #agregate decoded_cell
        rule:
            input: expand("results/"+rep+"/decoded/barcodes_unfiltered_no_neg_ctrl_ch_{{ch}}_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells),
            wildcard_constraints: cell="\d+", pos="\d+", dr="\d"
            output: "results/{rep}/decoded/barcodes_unfiltered_ch_{{ch}}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos,rep=rep)
            resources: mem_mb = 5000
            script: "../scripts/aggregate_sum_stats.py"

        rule:
            input: expand("results/"+rep+"/decoded/points_no_neg_ctrl_ch_{{ch}}_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells),
            wildcard_constraints: cell="\d+", pos="\d+", dr="\d"
            output: "results/{rep}/decoded/points_no_neg_ctrl_ch_{{ch}}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos,rep=rep)
            resources: mem_mb = 5000
            script: "../scripts/aggregate_sum_stats.py"
        
        rule:
            input: expand("results/"+rep+"/decoded/barcodes_unfiltered_w_neg_ctrl_ch_{{ch}}_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells),
            wildcard_constraints: cell="\d+", pos="\d+", dr="\d"
            output: "results/{rep}/decoded/barcodes_unfiltered_w_neg_ctrl_ch_{{ch}}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos,rep=rep)
            resources: mem_mb = 5000
            script: "../scripts/aggregate_sum_stats.py"

        rule:
            input: expand("results/"+rep+"/sum_stats/sum_stats_unfiltered_ch_{{ch}}_pos_"+str(pos)+"_cell{cell}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells),

            output: "results/{rep}/sum_stats/sum_stats_unfiltered_ch_{{ch}}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos,rep=rep)
            wildcard_constraints: pos="\d+"
            resources: mem_mb = 5000
            script: "../scripts/aggregate_sum_stats.py"



rule find_codepaths:
    #input: "aligned_dots_swp_h16_h61.csv", "codebooks/E2019_cb_all_control.txt"
    input: "results/rep{rep}/cell_dots/cell_dots_ch_{ch}_pos_{pos}_cell{cell}.csv",
           "results/full_codebooks/codebook_ch_{ch}.csv",
           "resources/codebooks/H.txt"
    #params: lf=min(lfs), wf = min(wfs), sf = min(sfs), dr = max(drops)
    params: lf = min(config['lfs']),
            wf = min(config['wfs']),
            sf = min(config['sfs']),
            dr = max(config['drops']),
	    rxy=config["rxy_ro"]
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    output: "results/rep{rep}/codepaths/cps_ch_{ch}_pos_{pos}_cell{cell}.csv"
    resources: mem_mb=20000
    script: "../scripts/find_save_codepath_candidates.jl"


rule choose_optimal_codepaths_w_neg_ctrl:
    input:
        "results/rep{rep}/codepaths/cps_ch_{ch}_pos_{pos}_cell{cell}.csv",
        "results/rep{rep}/cell_dots/cell_dots_ch_{ch}_pos_{pos}_cell{cell}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv",
        "resources/codebooks/H.txt"
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config['rxy_ro']
    group: "decode"
    output:
        #"results/points_lf{lf}_wf{wf}_sf{sf}_mw{dr}.csv",
        "results/rep{rep}/decoded/barcodes_unfiltered_w_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/sum_stats/sum_stats_w_neg_ctrl_unfiltered_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/decoded/points_w_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/dense_discarded_cpaths/dense_cpaths_w_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 40000
    script: "../scripts/choose_cpaths_from_saved_candidates.jl"


rule choose_optimal_codepaths_no_neg_ctrl:
    input:
        "results/rep{rep}/codepaths/cps_ch_{ch}_pos_{pos}_cell{cell}.csv",
        "results/rep{rep}/cell_dots/cell_dots_ch_{ch}_pos_{pos}_cell{cell}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv",
        "resources/codebooks/H.txt"
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config["rxy_ro"]#, rz = config["rz_ro"]
    output:
        "results/rep{rep}/decoded/barcodes_unfiltered_no_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/sum_stats/sum_stats_unfiltered_no_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/decoded/points_no_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/dense_discarded_cpaths/dense_cpaths_no_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", cell="\d+"
    resources: mem_mb = 12000
    script: "../scripts/choose_cpaths_from_saved_candidates_no_neg_ctrl.jl"

rule reconcile_sum_stats:
    input:
        "results/rep{rep}/sum_stats/sum_stats_w_neg_ctrl_unfiltered_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/sum_stats/sum_stats_unfiltered_no_neg_ctrl_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv"
    output: "results/rep{rep}/sum_stats/sum_stats_unfiltered_ch_{ch}_pos_{pos}_cell{cell}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    resources: mem_mb=1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/estimate_fdr.py"


rule filter_pos_ch_barcodes:
    input:
        "results/rep{rep}/decoded/barcodes_unfiltered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/full_codebooks/codebook_ch_{ch}.csv"
    params: filt_size=config["min_filter_size"], filt_prop=config["filter_prop_neg_control"]
    group: "decode"
    output:
        "results/rep{rep}/decoded/barcodes_filtered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/rep{rep}/sum_stats/sum_stats_filtered_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb=1000
    script: "../scripts/filter_barcodes.py"

rule get_fov_cnts:
    input: "results/rep{rep}/decoded/barcodes_{filtered}_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    output: "results/rep{rep}/decoded/counts_{filtered}_ch_{ch}_pos_{pos}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb=1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/get_fov_cnts.py"

rule aggregate_pos_ch_sum_stats:
    input: expand("results/rep{{rep}}/sum_stats/sum_stats_{{filtered}}_ch_{{ch}}_pos_{{pos}}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv",lf=config['lfs'],wf=config['wfs'],sf=config['sfs'],dr=config['drops'])
    output: "results/rep{rep}/pos_ch_summaries/decode_summary_stats_{filtered}_ch_{ch}_pos_{pos}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/aggregate_sum_stats.py"



# define rule for aggregating cell barcode counts and summary statistics across positions in for loop to account for differening number of positions in replicates
for rep in config['positions']:
    #aggregate_fov_counts
    rule:
        input: expand("results/"+rep+"/decoded/counts_{{filtered}}_ch_{ch}_pos_{pos}_lf{{lf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", pos=config['positions'][rep], ch=config['channels'])
        output: "results/"+rep+"/cell_barcode_counts/cell_barcode_counts_{filtered}_lf{lf}_wf{wf}_sf{sf}_dr{dr}.csv"
        resources: mem_mb=1000
        conda: "../envs/pdnp.yaml"
        script: "../scripts/aggregate_sum_stats.py"


    #rule aggregate_sum_stats:
    rule:
        input: expand("results/"+rep+"/pos_ch_summaries/decode_summary_stats_{{filtered}}_ch_{ch}_pos_{pos}.csv",ch=config['channels'],pos=config['positions'][rep],sf=config['sfs'])
        output: "results/"+rep+"/decode_summary_stats_{filtered}.csv"
        resources: mem_mb = 1000
        conda: "../envs/pdnp.yaml"
        script: "../scripts/aggregate_sum_stats.py"

rule all_reps_sum_stats:
    input: expand("results/{rep}/decode_summary_stats_{{filtered}}.csv",rep=config['positions']) #config['positions'] has a lists of the reps, then positions under the reps
    output: "results/decode_summary_stats_{filtered}.csv"
    resources: mem_mb=1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/get_all_rep_sum_stats.py"
