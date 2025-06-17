rule fit_beads:
    input: "results/rep{rep}/ch_pngs/{img_type}_ch_{ch}_pos_{pos}.png"
    params: sigma_lb = config['fit_beads']['sigma_lb'],
            sigma_ub = config['fit_beads']['sigma_ub'],
            final_loss_improvement = config['fit_beads']['final_loss_improvement'],
            min_weight = config['fit_beads']['min_weight'],
            max_iters = config['max_iters'],
            max_cd_iters = config['max_cd_iters'],
            min_allowed_separation = config['min_allowed_separation']
    output: "results/rep{rep}/fit_beads/{img_type}_ch_{ch}_pos_{pos}.csv"
    group: "fit_beads"
    resources: mem_mb = 4000
    script: "../scripts/tile_find_beads.jl"

rule aggregate_bead_fits:
    input: expand("results/rep{{rep}}/fit_beads/HybCycle_{hyb}_ch_{{ch}}_pos_{{pos}}.csv", hyb = hybs)
    output: "results/rep{rep}/fit_beads/all_hyb_bright_dots_ch_{ch}_pos_{pos}.csv"
    resources: mem_mb = 2000
    group: "align"
    conda: "../envs/pdnp.yaml"
    script: "../scripts/aggregate_dots.py"
