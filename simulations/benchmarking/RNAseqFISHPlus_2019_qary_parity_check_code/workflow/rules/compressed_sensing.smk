for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1
    pc_range = list(range(1, code_q + 1)) if code[:7] == "seqFISH" else list(range(1, code_q))

    rule:
        name: f"cs_get_candidate_coords_{code}"
        input:
            "resources/codebooks/" + code + "cb.csv",
            "workflow/envs/cs_julia_environment/",
            "results/cs_julia_environment_installed.txt",
            expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=pc_range
            )
        output:
            "results/compressed_sensing/cand_dots/cdots_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
        params:
            fov_width=config["roi_width"] + 2 * config["roi_pad"],
            grid_size=config["grid_size"],
            min_intensity=config["min_weight"],
            code=code
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"])
        benchmark: "results/benchmarks/cs_spot_detection/" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
        script: "../scripts/get_candidate_coords.jl"


rule cs_find_cand_cpaths:
    input:
        "resources/codebooks/{code}cb.csv",
        "resources/codebooks/{code}_H.csv",
        "results/compressed_sensing/cand_dots/cdots_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        "workflow/envs/cs_julia_environment/",
        "results/cs_julia_environment_installed.txt"
    output:
        "results/compressed_sensing/cand_cpaths/cand_cpaths_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        sr=config["cs_search_radius"],
        grid_size=config["grid_size"]
    wildcard_constraints:
        rep=r"\d+",
        code="|".join(codes_to_simulate["code"])
    benchmark: "results/benchmarks/cs_cand_cpaths/{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
    script: "../scripts/get_cand_cpaths.jl"


rule cs_find_cand_cpaths_with_negctrls:
    input:
        "resources/codebooks/{code}cb.csv",
        "resources/codebooks/{code}_H.csv",
        "results/compressed_sensing/cand_dots/cdots_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        "workflow/envs/cs_julia_environment/",
        "results/cs_julia_environment_installed.txt"
    output:
        "results/compressed_sensing/cand_cpaths/nc_cand_cpaths_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        sr=config["cs_search_radius"],
        grid_size=config["grid_size"]
    wildcard_constraints:
        rep=r"\d+",
        code="|".join(codes_to_simulate["code"])
    benchmark: "results/benchmarks/cs_cand_cpaths_nc/{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
    script: "../scripts/get_cand_cpaths.jl"


import platform

for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1
    pc_range = list(range(1, code_q + 1)) if code[:7] == "seqFISH" else list(range(1, code_q))

    rule:
        name: f"cs_l0learn_fit_{code}"
        input:
            "results/compressed_sensing/cand_cpaths/nc_cand_cpaths_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "results/compressed_sensing/cand_dots/cdots_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "resources/codebooks/" + code + "cb.csv",
            expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=pc_range,
                code=code
            )
        output:
            "results/compressed_sensing/decoded_all_ranked/nc_decoded_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "results/compressed_sensing/lambda_stats/nc_lambda_stats_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "results/compressed_sensing/decoded/nc_decoded_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
            "results/compressed_sensing/decoded_perlambda/nc_perlambda_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
        params:
            code=code
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"])
        benchmark: "results/benchmarks/cs_l0learn/" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
        # No conda package provides L0Learn for win-64 (see r_env.yaml), so on Windows we
        # skip conda entirely and run against the system R install, which already has the
        # required packages (L0Learn, glmnet, Matrix, png, tiff).
        conda: None if platform.system() == "Windows" else "../envs/r_env.yaml"
        script: "../scripts/l0learn_barcode_hyperstack_fit_sparseX.R"
