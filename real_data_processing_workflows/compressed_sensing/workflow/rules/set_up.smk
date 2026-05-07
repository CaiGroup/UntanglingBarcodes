rule install_julia_environment:
    input: "workflow/envs/julia_environment/"
    output: "results/julia_environment_installed.txt"
    script: "../scripts/install_julia_environment.jl"