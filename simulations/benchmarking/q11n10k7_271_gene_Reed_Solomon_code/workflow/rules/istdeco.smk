import itertools

for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1

    # One job per condition sweeps every (Q, tau) combo itself: model.run() (the expensive
    # ISTDeco optimization) doesn't depend on either threshold, so it now runs once per
    # condition instead of once per (Q, tau) pair, and the cheap per-gene thresholding is
    # applied for every combo against that single result.
    istdeco_combos = list(itertools.product(config["istdeco_threshold_Q"], config["istdeco_threshold_percentile"]))

    rule:
        name: f"istdeco_decode_{code}"
        input:
            "resources/codebooks/" + code + "cb.csv",
            *expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=list(range(1, code_q))
            )
        output:
            [
                "results/istdeco/decoded/decoded_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_Q"
                + str(Q) + "_tau" + str(tau) + ".csv"
                for Q, tau in istdeco_combos
            ]
        params:
            sigma=config["sigma"],
            niter=config["istdeco_niter"],
            suppress_radius=config["istdeco_suppress_radius"],
            b=config["istdeco_b"],
            device=config["istdeco_device"],
            istdeco_path="workflow/scripts/istdeco_lib",
            combos=istdeco_combos
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"])
        #benchmark: "results/benchmarks/istdeco/" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
        conda: "../envs/istdeco.yaml"
        script: "../scripts/run_istdeco.py"
