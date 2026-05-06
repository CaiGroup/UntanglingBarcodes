rule fit_beads:
    input: 
        "resources/imgs/pos{pos}/HybCycle_{hyb}_pos_{pos}_z_5.png", #"resources/imgs/Cy3B_Pos{pos}_png/pos{pos}_hyb_{hyb}.png",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
    params: sigma_lb = config['fit_beads']['sigma_lb'],
            sigma_ub = config['fit_beads']['sigma_ub'],
            final_loss_improvement = config['fit_beads']['final_loss_improvement'],
            min_weight = config['fit_beads']['min_weight'],
            max_iters = config['max_iters'],
            max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation']
    output: "results/fit_beads/pos_{pos}_hyb_{hyb}.csv"
    group: "fit_beads"
    benchmark: "benchmarks/fit_beads/pos_{pos}_hyb_{hyb}.txt"
    resources: mem_mb = 4000
    script: "../scripts/tile_find_beads.jl"

rule fit_beads_init_reference:
    input: 
        "resources/imgs/pos{pos}/initial_background_MMStack_Pos{pos}_z_5.png",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
    params: sigma_lb = config['fit_beads']['sigma_lb'],
            sigma_ub = config['fit_beads']['sigma_ub'],
            final_loss_improvement = config['fit_beads']['final_loss_improvement'],
            min_weight = config['fit_beads']['min_weight'],
            max_iters = config['max_iters'],
            max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation']
    output: "results/fit_beads/beads_initial_pos{pos}.csv"
    group: "fit_beads"
    benchmark: "benchmarks/fit_beads/beads_initial_pos{pos}.txt"
    resources: mem_mb = 4000
    script: "../scripts/tile_find_beads.jl"

rule fit_beads_final_reference:
    input: 
        "resources/imgs/pos{pos}/final_background_MMStack_Pos{pos}_z_5.png",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt"
    params: sigma_lb = config['fit_beads']['sigma_lb'],
            sigma_ub = config['fit_beads']['sigma_ub'],
            final_loss_improvement = config['fit_beads']['final_loss_improvement'],
            min_weight = config['fit_beads']['min_weight'],
            max_iters = config['max_iters'],
            max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation']
    output: "results/fit_beads/beads_final_pos{pos}.csv"
    group: "fit_beads"
    benchmark: "benchmarks/fit_beads/beads_final_pos{pos}.txt"
    resources: mem_mb = 4000
    script: "../scripts/tile_find_beads.jl"

rule aggregate_bead_fits:
    input: expand("results/fit_beads/pos_{{pos}}_hyb_{hyb}.csv", hyb = range(100)) #range(1,101))
    output: "results/fit_beads/all_hyb_bright_dots_pos_{pos}.csv"
    resources: mem_mb = 2000
    group: "align"
    benchmark: "benchmarks/fit_beads/all_hyb_bright_dots_pos_{pos}.txt"
    conda: "../envs/pdnp.yaml"
    script: "../scripts/aggregate_dots.py"
