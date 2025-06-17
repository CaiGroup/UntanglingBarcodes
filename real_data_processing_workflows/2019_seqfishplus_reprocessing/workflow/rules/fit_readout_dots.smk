rule tile_fit_dots:
    input: "results/rep{rep}/ims_bg_sub/HybCycle_{hyb}_ch_{ch}_pos_{pos}.png",
           "results/rep{rep}/Shifted_Labeled_Images/hyb_{hyb}_ch_{ch}_pos_{pos}.png"
    params: sigma_lb = config['sigma_lb'],
            sigma_ub = config['sigma_ub'],
            final_loss_improvement = config['final_loss_improvement'],
            min_weight = config['min_weight'],
            max_iters = config['max_iters'],
            max_cd_iters = config['max_cd_iters']
    output: "results/rep{rep}/fit_dots/HybCycle_{hyb}_ch_{ch}_pos_{pos}_w_duplicates.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+", hyb="\d+"
    group: "fit_dots"
    resources: mem_mb = 3000
    script: "../scripts/tile_find_dots.jl"


rule remove_duplicates:
    #input: "ims_bg_sub/HybCycle_{hyb}_ch_{ch}.png", "fit_dots/HybCycle_{hyb}_mw_{mw}_ch_{ch}_w_duplicates.csv"
    input:
        "results/rep{rep}/ims_bg_sub/HybCycle_{hyb}_ch_{ch}_pos_{pos}.png",
        "results/rep{rep}/fit_dots/HybCycle_{hyb}_ch_{ch}_pos_{pos}_w_duplicates.csv",
        "results/rep{rep}/Shifted_Labeled_Images/hyb_{hyb}_ch_{ch}_pos_{pos}.png"
    #params: sigma_lb = 1.0*488/647, sigma_ub = 1.4*488/647, min_allowed_separation = 2.0
    params:
            sigma_lb = config['sigma_lb'],
            sigma_ub = config['sigma_ub'],
            min_allowed_separation = config['min_allowed_separation']
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    group: "fit_dots"
    output: "results/rep{rep}/fit_dots/HybCycle_{hyb}_ch_{ch}_pos_{pos}.csv"
    resources: mem_mb = 2000
    script: "../scripts/remove_duplicates.jl"


rule aggregate_hyb_points:
    #input: expand("fit_dots/{hyb}.csv", hyb = hybs)
    input: expand("results/rep{{rep}}/fit_dots/HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}.csv", hyb = hybs)
    output: "results/rep{rep}/fit_dots/all_hyb_dots_ch_{ch}_pos_{pos}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 2000
    conda: "../envs/pdnp.yaml"
    group: "decode"
    script: "../scripts/aggregate_dots.py"


rule apply_offsets:
    input: "results/rep{rep}/alignment/offsets_ch_{ch}_pos_{pos}.csv",
           "results/rep{rep}/fit_dots/all_hyb_dots_ch_{ch}_pos_{pos}.csv"
    output: "results/rep{rep}/alignment/aligned_dots_ch_{ch}_pos_{pos}.csv"
    wildcard_constraints: rep="\d+", ch="\d+", pos="\d+"
    resources: mem_mb = 2000
    group: "fit_dots"
    conda: "../envs/pdnp.yaml"
    script: "../scripts/apply_alignments.py"
