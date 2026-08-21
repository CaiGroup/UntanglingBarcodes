rule install_cs_julia_environment:
    input: "workflow/envs/cs_julia_environment"
    output: "results/cs_julia_environment_installed.txt"
    resources: mem_mb=10000
    script: "../scripts/install_julia_environment.jl"


rule generate_ground_truth:
    input:
        "resources/codebooks/{code}cb.csv",
        "workflow/envs/cs_julia_environment/",
        "results/cs_julia_environment_installed.txt"
    params:
        roi_width=config["roi_width"],
        roi_pad=config["roi_pad"],
        sigma=config["sigma"],
        rstdv=lambda wildcards: float(wildcards.rstdv),
        pbind_primary=config["pbind_primary"]
    output:
        "results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    wildcard_constraints:
        rep=r"\d+",
        code="|".join(codes_to_simulate["code"]),
        rstdv=r"[\d.]+"
    script: "../scripts/gen_truth.jl"


for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1

    rule:
        name: f"generate_sim_data_{code}"
        input:
            "resources/codebooks/" + code + "cb.csv",
            "results/truth/ground_truth_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "workflow/envs/cs_julia_environment/",
            "results/cs_julia_environment_installed.txt"
        params:
            fov_width=config["roi_width"] + 2 * config["roi_pad"],
            sigma=config["sigma"],
            rstdv=lambda wildcards: float(wildcards.rstdv),
            dot_intensity=config["dot_intensity"],
            pbind_primary=config["pbind_primary"],
            pbind_secondary=config["pbind_readout"]
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"]),
            rstdv=r"[\d.]+"
        output:
            "results/truth/ground_truth_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
            "results/dots/dots_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=list(range(1, code_q))
            )
        script: "../scripts/gen_sim_data.jl"
