rule generate_full_cb:
    input: 
        "resources/codebooks/channel_{ch}.csv",
        "resources/codebooks/H.txt",
        "workflow/envs/julia_environment",
        "results/julia_environment_installed.txt" 
    params: n_neg_cntrl_cws = config["n_neg_cntrl_cws"]
    output: "results/full_codebooks/codebook_ch_{ch}.csv"
    resources: mem_mb = 100
    group: "gen_full_cb"
    script: "../scripts/gen_cntrl_cb.jl"
