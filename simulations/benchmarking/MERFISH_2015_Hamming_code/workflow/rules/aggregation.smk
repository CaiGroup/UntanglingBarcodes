methods_other = ["bardensr", "istdeco", "compressed_sensing", "spotsfirst"]

rule aggregate_all_for_code:
    input:
        expand(
            "results/scores/{method}/score_{{code}}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            method=methods_other,
            rep=config["reps"],
            nbarcodes=config["nbarcodes"],
            nnonspec=config["nnonspec"],
            pdrop=config["pdrop"],
            rstdv=config["rstdv"]
        )
    output:
        "results/all_methods_{code}.csv"
    conda: "../envs/python_base.yaml"
    script: "../scripts/aggregate_all_methods.py"


rule aggregate_benchmarks:
    input:
        expand(
            "results/all_methods_{code}.csv",
            code=codes_to_simulate["code"]
        )
    output:
        "results/benchmarks/all_benchmarks.csv"
    conda: "../envs/python_base.yaml"
    script: "../scripts/aggregate_benchmarks.py"


rule plot_computational_cost:
    input:
        "results/benchmarks/all_benchmarks.csv"
    output:
        "results/plots/computational_cost.png",
        "results/plots/computational_cost.pdf"
    conda: "../envs/python_base.yaml"
    script: "../scripts/plot_computational_cost.py"


rule plot_sens_fdr:
    input:
        "results/all_methods_{code}.csv"
    output:
        "results/plots/sens_fdr_vs_density_{code}.png",
        "results/plots/sens_fdr_vs_density_{code}.pdf"
    conda: "../envs/python_base.yaml"
    script: "../scripts/plot_sens_fdr_vs_density.py"
