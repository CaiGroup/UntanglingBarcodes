rule find_alignment_parameters:
    input: "results/fit_beads/beads_initial_pos{pos}.csv",
           "results/fit_beads/beads_final_pos{pos}.csv",
           "results/fit_beads/all_hyb_bright_dots_pos_{pos}.csv"
    output: "results/alignment/params_pos_{pos}.csv"
    benchmark: "benchmarks/alignment/params_pos_{pos}.txt"
    conda: "../envs/align.yaml"
    resources: mem_mb=500
    group: "align"
    script: "../scripts/find_alignment_params.py"

rule find_alignments:
    input: "results/fit_beads/beads_initial_pos{pos}.csv",
           "results/fit_beads/beads_final_pos{pos}.csv",
           "results/fit_beads/all_hyb_bright_dots_pos_{pos}.csv",
           "results/alignment/params_pos_{pos}.csv"
    output: "results/alignment/offsets_pos_{pos}.csv",
            "results/alignment/matches_pos_{pos}.csv",
            "results/alignment/loov_errors_pos_{pos}.csv"#,"alignment/pnts_no_fm.csv"
    benchmark: "benchmarks/alignment/offsets_pos_{pos}.txt"
    params:
        min_fm_hyb_matches = config['min_fm_hyb_matches'],
        outlier_sd_thresh = config['outlier_sd_thresh'],
        max_lat_offset = config['max_lat_offset'],
        set_xy_search_error = config['set_xy_search_error'] 
    conda: "../envs/align.yaml"
    resources: mem_mb=500
    #group: "align"
    script: "../scripts/align.py"

rule apply_alignments:
    input: "results/alignment/offsets_pos_{pos}.csv",
           "results/fit_dots/all_ro_pos_{pos}_z_{z}.csv", #"results/cand_dots/cand_dots_pos_{pos}_cell_{cell}_mpcd{mpcd}.csv",
           "resources/codebooks/block_pseuedocolor_hyb_table.csv",
           labeled_image = "results/filtered_labeled_imgs/lbld_img_pos{pos}_z_{z}.png"
    group: "find_dots"
    output: "results/aligned_dots/aligned_dots_pos_{pos}_z_{z}.csv"
    benchmark: "benchmarks/filtered_labeled_imgs/lbld_img_pos{pos}_z_{z}.txt"
    script: "../scripts/apply_alignments.py"