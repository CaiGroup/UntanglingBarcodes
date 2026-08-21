for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1
    pc_range = list(range(1, code_q + 1)) if code[:7] == "seqFISH" else list(range(1, code_q))

    rule:
        name: f"graph_iss_decode_{code}"
        input:
            "resources/codebooks/" + code + "cb.csv",
            *expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=pc_range
            )
        output:
            "results/graph_iss/decoded/decoded_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_k1{k1}.csv"
        params:
            spot_min_distance=config["graph_iss_spot_min_distance"],
            spot_threshold_rel=config["graph_iss_spot_threshold_rel"],
            spatial_radius=config["graph_iss_spatial_radius"],
            k1=lambda wildcards: float(wildcards.k1),
            code=code
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"]),
            k1=r"[\d.]+"
        #benchmark: "results/benchmarks/graph_iss/" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
        conda: "../envs/graph_iss.yaml"
        script: "../scripts/run_graph_iss.py"
