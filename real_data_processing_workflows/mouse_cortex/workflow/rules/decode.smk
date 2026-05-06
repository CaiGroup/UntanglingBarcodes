import skimage
import numpy as np

for pos in config["positions"]:
    #rule split_cells
    lbld_img = skimage.io.imread("resources/labeled_stacks/segmentation_MMStack_Pos{pos}.tif".format(pos=pos))
    np.unique(lbld_img)
    cells = [cellID for cellID in np.unique(lbld_img) if cellID != 0]
    rule:
        input: "results/aligned_dots/aligned_dots_pos_{pos}_z_{{z}}.csv".format(pos=pos)
        output: expand("results/cell_dots/cell_dots_pos"+str(pos)+"_z_{{z}}_cell{cell}.csv", cell=cells)
        #benchmark: "benchmarks/cell_dots/cell_dots_pos"+str(pos)+".txt"
        wildcard_constraints: cell="\d+", pos="\d+", dr="\d", z="\d+"
        resources: mem_mb = 5000
        script: "../scripts/split_cells.py"

    #agregate decoded_cell
    rule:
        input: expand("results/decoded/barcodes_unfiltered_pos_"+str(pos)+"_z_{{z}}_cell{cell}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells), 
        wildcard_constraints: cell="\d+", pos="\d+", dr="\d", z="\d+"
        output: "results/decoded/barcodes_unfiltered_pos_{pos}_z_{{z}}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos)
        benchmark: "benchmark/decoded/barcodes_unfiltered_pos_{pos}_z_{{z}}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos)
        resources: mem_mb = 5000
        script: "../scripts/aggregate_sum_stats.py"

    rule:
        input: expand("results/sum_stats/sum_stats_unfiltered_pos_"+str(pos)+"_z_{{z}}_cell{cell}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv", cell=cells), 
        output: "results/sum_stats/sum_stats_unfiltered_pos_{pos}_z_{{z}}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv".format(pos=pos)
        benchmark: f"benchmarks/sum_stats/sum_stats_unfiltered_pos_{pos}_z_{{z}}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv.txt"
        wildcard_constraints: pos="\d+", z="\d+"
        resources: mem_mb = 5000
        script: "../scripts/aggregate_sum_stats.py"


rule find_codepaths:
    input: 
           "results/cell_dots/cell_dots_pos{pos}_z_{z}_cell{cell}.csv",
           "resources/codebooks/full_RS_q11_k7_cb.csv",
           "resources/codebooks/RS_q11_k7_H.csv",
           "workflow/envs/julia_environment",
           "results/julia_environment_installed.txt"

    #params: lf=min(lfs), wf = min(wfs), sf = min(sfs), dr = max(drops)
    params: lf = min(config['lfs']),
            zf = min(config['zfs']),
            wf = min(config['wfs']),
            sf = min(config['sfs']),
            dr = max(config['drops']),
            rxy = config["rxy_ro"], 
            rz = config["rz_ro"]
    wildcard_constraints: ch="\d+", pos="\d+", z= "\d+"
    output: "results/codepaths/cps_pos_{pos}_z_{z}_cell{cell}.csv"
    benchmark: "benchmarks/codepaths/cps_pos_{pos}_z_{z}_cell{cell}.txt"
    group: "find codepaths"
    resources: mem_mb=20000
    script: "../scripts/find_save_codepath_candidates.jl"


rule choose_optimal_codepaths_w_neg_ctrl:
    input:
        "results/codepaths/cps_pos_{pos}_z_{z}_cell{cell}.csv",
        "results/cell_dots/cell_dots_pos{pos}_z_{z}_cell{cell}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv",
        "resources/codebooks/RS_q11_k7_H.csv",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config["rxy_ro"], rz = config["rz_ro"]
    group: "decode"
    output:
        "results/decoded/barcodes_unfiltered_pos_{pos}_z_{z}cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_w_neg_ctrl_unfiltered_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/decoded/points_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/dense_discarded_cpaths/dense_cpaths__pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    benchmark: "benchmarks/decoded/barcodes_unfiltered_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.txt"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", cell="\d+", z="\d+"
    resources: mem_mb = 40000
    script: "../scripts/choose_cpaths_from_saved_candidates.jl"

rule choose_optimal_codepaths_no_neg_ctrl:
    input:
        "results/codepaths/cps_pos_{pos}_z_{z}_cell{cell}.csv",
        "results/cell_dots/cell_dots_pos{pos}_z_{z}_cell{cell}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv",
        "resources/codebooks/RS_q11_k7_H.csv",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
    params: skip_thresh = config["skip_thresh"], skip_density_thresh=config["skip_density_thresh"], rxy = config["rxy_ro"], rz = config["rz_ro"]
    output:
        "results/decoded/barcodes_unfiltered_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_unfiltered_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/decoded/points_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf_{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/dense_discarded_cpaths/dense_cpaths_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    benchmark: "benchmarks/decoded/barcodes_unfiltered_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.txt"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", cell="\d+", z="\d+"
    resources: mem_mb = 20000
    script: "../scripts/choose_cpaths_from_saved_candidates_no_neg_ctrl.jl"

rule reconcile_sum_stats:
    input:
        "results/sum_stats/sum_stats_w_neg_ctrl_unfiltered_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "results/sum_stats/sum_stats_unfiltered_no_neg_ctrl_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv",
        "resources/codebooks/full_RS_q11_k7_cb.csv"
    output: "results/sum_stats/sum_stats_unfiltered_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    benchmark: "benchmarks/sum_stats/sum_stats_unfiltered_pos_{pos}_z_{z}_cell{cell}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.txt"
    resources: mem_mb=1000
    conda: "../envs/pdnp.yaml"
    script: "../scripts/estimate_fdr.py"


rule aggregate_sum_stats:
    input: expand("results/sum_stats/sum_stats_unfiltered_pos_{{pos}}_z_{z}_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv", lf=config["lfs"], zf=config["zfs"], wf=config["wfs"], sf=config["sfs"], dr=config["drops"], pos=config["positions"], z=config['zslices'])
    output: "results/decode_summary_stats_pos_{pos}.csv"
    benchmark: "benchmarks/decode_summary_stats_pos_{pos}.txt"
    script: "../scripts/aggregate_sum_stats_z.py"


rule aggregate_decoded_position:
    input:
        "resources/pos.pos",
        expand('results/decoded/barcodes_unfiltered_pos_{pos}_lf{{lf}}_zf{{zf}}_wf{{wf}}_sf{{sf}}_dr{{dr}}.csv', pos=config['positions'])
    output:"results/decoded_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.csv"
    benchmark: "benchmarks/decoded_lf{lf}_zf{zf}_wf{wf}_sf{sf}_dr{dr}.txt"
    script: "../scripts/aggregate_positions.py"