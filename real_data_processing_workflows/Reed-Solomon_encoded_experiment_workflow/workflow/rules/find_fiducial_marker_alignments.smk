rule find_alignment_parameters:
    input: "results/fm_fits/HybCycle_0_ch{ch}_pos{pos}.csv",
           "results/fm_fits/all_fms_pos_{pos}.csv"
    output: "results/alignment/params_pos_{pos}.csv"
    conda: "../envs/pdnp.yaml"
    resources: mem_mb=500
    group: "align"
    script: "../scripts/find_alignment_params.py"

rule find_alignments:
    input: "results/fm_fits/chromatic_aberration_ch{ch}_pos{pos}.csv",
           "results/fm_fits/all_fms_pos_{pos}.csv"

    output: "results/alignment/offsets_ch{ch}_pos_{pos}.csv",
            "results/alignment/matches_ch{ch}_pos_{pos}.csv",
            "results/alignment/loov_errors_ch{ch}_pos_{pos}.csv"#,"alignment/pnts_no_fm.csv"
    conda: "../envs/align.yaml"
    resources: mem_mb=500
    group: "align"
    script: "../scripts/align.py"

rule match_cross_channel_fiducials:
    input: expand("results/fm_fits/chromatic_aberration_ch{ch}_pos{{pos}}.csv", ch=config["channels"]),
    output: "results/alignment/cross_channel_matched_reference_dots_pos{pos}.csv",
    conda: "../envs/align.yaml"
    resources: mem_mb=500
    script: "../scripts/match_cross_channel_fiducials.py"

rule register_to_reference:
    input:
       "results/alignment/cross_channel_matched_reference_dots_pos{pos}.csv",
       "results/fit_dots/all_ro_pos_{pos}.csv",
       "results/alignment/matches_ch2_pos_{pos}.csv",
       "results/alignment/matches_ch1_pos_{pos}.csv",
       "results/alignment/matches_ch0_pos_{pos}.csv",
       "results/alignment/offsets_ch2_pos_{pos}.csv",
       "results/alignment/offsets_ch1_pos_{pos}.csv",
       "results/alignment/offsets_ch0_pos_{pos}.csv",
       "resources/codebooks/block_pseudocolor_table.csv",
       "results/filtered_labeled_imgs/lbld_img_pos{pos}.png" #"results/labeled_stacks/labeled_stack_pos{pos}.tif"
    output: "results/alignment/ref_cross_channel_registered_dots_pos{pos}.csv", "results/alignment/ref_cross_channel_loocv_pos{pos}.csv"
    conda: "../envs/register.yaml"
    resources: mem_mb=5000
    script: "../scripts/register_affine_to_ref.py"


"""
rule register_to_reference:
    input:
       "results/alignment/cross_channel_matched_reference_dots_pos{pos}.csv",
       "results/fit_dots/all_ro_pos_{pos}.csv",
       "results/alignment/matches_ch488_pos_{pos}.csv",
       "results/alignment/matches_ch561_pos_{pos}.csv",
       "results/alignment/matches_ch640_pos_{pos}.csv",
       "results/alignment/offsets_ch488_pos_{pos}.csv",
       "results/alignment/offsets_ch561_pos_{pos}.csv",
       "results/alignment/offsets_ch640_pos_{pos}.csv",
       "resources/codebooks/block_pseudocolor_table.csv",
       "results/labeled_stacks/labeled_stack_pos{pos}.tif"
    output: "results/alignment/ref_cross_channel_registered_dots_pos{pos}.csv"
    resources: mem_mb=500
    script: "../scripts/register_affine_to_ref.py"
"""
