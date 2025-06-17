
rule find_alignment_parameters:
    input: "results/rep{rep}/fit_beads/beads_initial_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/fit_beads/beads_final_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/fit_beads/all_hyb_bright_dots_ch_{ch}_pos_{pos}.csv"
    output: "results/rep{rep}/alignment/params_ch_{ch}_pos_{pos}.csv"
    resources: mem_mb = 2000
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    group: "align"
    conda: '../envs/align.yaml'
    script: "../scripts/find_alignment_params.py"

rule find_alignments:
    input: "results/rep{rep}/fit_beads/beads_initial_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/fit_beads/beads_final_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/fit_beads/all_hyb_bright_dots_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/alignment/params_ch_{ch}_pos_{pos}.csv"
    params:
        min_fm_hyb_matches = config['min_fm_hyb_matches'],
        outlier_sd_thresh = config['outlier_sd_thresh'],
        max_lat_offset = config['max_lat_offset'],
        set_xy_search_error = config['set_xy_search_error'] #config['r_xy'], config['r_z']
    output: "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
            "results/rep{rep}/alignment/matches_ch_{ch}_pos_{pos}.csv",
            "results/rep{rep}/alignment/loov_errors_ch_{ch}_pos_{pos}.csv"#,"alignment/pnts_no_fm.csv"
    resources: mem_mb = 2000
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    group: "align"
    conda: '../envs/align.yaml'
    script: "../scripts/align.py"
