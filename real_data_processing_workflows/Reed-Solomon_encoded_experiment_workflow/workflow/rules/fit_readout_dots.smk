
rule fit_ro:
    input: 
        "results/ims_bg_sub/HybCycle_{hyb}_ch{ch}_pos_{pos}.png",
        "results/filtered_labeled_imgs/lbld_img_pos{pos}.png"
        #"results/labeled_stacks/labeled_stack_pos{pos}.tif"
    #input: "results/ch_stacks/HybCycle_{hyb}_ch_{ch}_pos{pos}.tif",
    params: sigma_xy_lb = config["sigma_xy_lb"], sigma_xy_ub = config["sigma_xy_ub"], sigma_z_lb = config["sigma_z_lb"], sigma_z_ub = config["sigma_z_ub"],
            final_loss_improvement = config["final_loss_improvement"],
            min_weight = config['min_weight'],
            max_iters = config['max_iters'], max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation'],
            tile_main_width = config["tile_main_width"], tile_overlap = config["tile_overlap"],
            tile_depth = config["tile_depth"], tile_depth_overhang = config["tile_depth_overhang"]
    wildcard_constraints: hyb="\d+", ch="\d+", pos="\d+"
    output:
        "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.csv"#,
        #"results/fit_dots/records_w_dups_HybCycle_{hyb}_ch_{ch}_pos{pos}.txt"
    resources: mem_mb=9000
    group: "fit_dots"
    script: "../scripts/tile_fit_3d.jl"


"""

rule fit_irafstarfinder:
    input: 
        "results/ims_bg_sub/HybCycle_{hyb}_ch{ch}_pos_{pos}.tif",
    output:
        "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.csv"
    script:
        "../scripts/fit_irafstarfinder.py"

rule fit_daostorm:
    input: 
        "results/ims_bg_sub/HybCycle_{hyb}_ch{ch}_pos_{pos}.tif",
        "resources/dao_params.xml"
    output:
        "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.hdf5"
    script:
        "../scripts/fit_daostorm.py"

rule daostorm_out_2_csv:
    input: "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.hdf5"
    output: "results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{pos}.csv"
    script: "../scripts/storm_output_2_csv.py"
"""

"""
rule get_intermediate_mw_fit_res:
    input: "results/fit_dots/records_w_dups_HybCycle_{hyb}_ch_{ch}_pos{pos}.csv"
    params: sigma_xy_lb = 1.0, sigma_xy_ub = 2.2, min_allowed_separation = config['min_allowed_separation'], mws = config["min_weights"]
    output: expand("results/fit_dots/HybCycle_{{hyb}}_ch_{{ch}}_pos{{pos}}_mw_{mw}.txt", mw=config["min_weights"])
    group: "fit_dots"
    script: "../scripts/get_intermediate_mw_res.jl"
"""

rule aggregate_ro_fits:
    #input: expand("fit_dots/{hyb}.csv", hyb = hybs)
    #input: expand("results/fit_dots/HybCycle_{hyb}_ch_{{ch}}_pos{{pos}}_mw_{{mw}}.txt", hyb = hybs)
    input: expand("results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{{pos}}.csv", hyb = hybs, ch=config["channels"])
    #output: "results/fit_dots/all_ro_ch_{ch}_pos_{pos}_mw_{mw}.csv"
    output: "results/fit_dots/all_ro_pos_{pos}.csv"
    wildcard_constraints: pos="\d+"
    conda: "../envs/pdnp.yaml"
    resources: mem_mb=2000
    script: "../scripts/aggregate_dots.py"

rule aggregate_smFISH_fits:
    input:
        "resources/codebooks/smFISH_codebook.csv",
         "results/filtered_labeled_imgs/lbld_img_pos{pos}.png", #"results/labeled_stacks/labeled_stack_pos{pos}.tif",
        expand("results/fit_dots/HybCycle_{hyb}_ch_{ch}_pos{{pos}}.csv", hyb = smFISH_hybs, ch=config["channels"])
    output: "results/fit_dots/all_smFISH_pos_{pos}.csv"
    conda: "../envs/pdnp.yaml"
    resources: mem_mb=2000
    script: "../scripts/aggregate_smFISH_fits.py"
