
rule fit_ro:
    input: 
        "results/ims_bg_sub/HybCycle_{hyb}_pos_{pos}_z_{z}.png",
        "results/filtered_labeled_imgs/lbld_img_pos{pos}_z_{z}.png",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
        #"results/labeled_stacks/labeled_stack_pos{pos}.tif"
    #input: "results/ch_stacks/HybCycle_{hyb}_ch_{ch}_pos{pos}.tif",
    params: sigma_xy_lb = config["sigma_xy_lb"], sigma_xy_ub = config["sigma_xy_ub"], sigma_z_lb = config["sigma_z_lb"], sigma_z_ub = config["sigma_z_ub"],
            final_loss_improvement = config["final_loss_improvement"],
            min_weight = config['min_weight'],
            max_iters = config['max_iters'], max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation'],
            tile_main_width = config["tile_main_width"], tile_overlap = config["tile_overlap"],
            tile_depth = config["tile_depth"], tile_depth_overhang = config["tile_depth_overhang"]
    wildcard_constraints: hyb="\d+", pos="\d+"
    output: "results/fit_dots/HybCycle_{hyb}_pos{pos}_z_{z}.csv"
    benchmark: "benchmarks/fit_dots/HybCycle_{hyb}_pos{pos}_z_{z}.txt"
    resources: mem_mb=9000
    group: "fit_dots"
    script: "../scripts/tile_fit.jl"


rule aggregate_ro_fits:
    #input: expand("fit_dots/{hyb}.csv", hyb = hybs)
    #input: expand("results/fit_dots/HybCycle_{hyb}_ch_{{ch}}_pos{{pos}}_mw_{{mw}}.txt", hyb = hybs)
    input: expand("results/fit_dots/HybCycle_{hyb}_pos{{pos}}_z_{{z}}.csv", hyb = list(range(100))) #range(1,101)))
    #output: "results/fit_dots/all_ro_ch_{ch}_pos_{pos}_mw_{mw}.csv"
    output: "results/fit_dots/all_ro_pos_{pos}_z_{z}.csv"
    benchmark: "benchmarks/fit_dots/all_ro_pos_{pos}_z_{z}.txt"
    wildcard_constraints: pos="\d+"
    conda: "../envs/pdnp.yaml"
    resources: mem_mb=2000
    script: "../scripts/aggregate_ro_dots.py"