import itertools

for code in codes_to_simulate["code"]:
    path = "resources/codebooks/{code}cb.csv".format(code=code)
    code_cb = pd.read_csv(path)
    code_q = len(np.unique(code_cb.iloc[:, 1:]))
    code_n = np.shape(code_cb)[1] - 1
    pc_range = list(range(1, code_q + 1)) if code[:7] == "seqFISH" else list(range(1, code_q))

    # One job per condition sweeps every (l1, peak) combo itself: estimate_density_iterative
    # (the expensive step) only depends on l1, so this reuses one density map across all peak
    # thresholds for that l1 instead of recomputing it once per (l1, peak) pair, and it also
    # collapses what was previously len(l1)*len(peak) separate jobs (each paying its own conda
    # activation + image load + TF startup) into a single job per condition.
    bardensr_combos = list(itertools.product(config["bardensr_l1_penalty"], config["bardensr_peak_thresh"]))

    rule:
        name: f"bardensr_decode_{code}"
        input:
            "resources/codebooks/" + code + "cb.csv",
            *expand(
                "results/images/hyb_im_r_{r}_pc_{pc}_" + code + "_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}.png",
                r=list(range(1, code_n + 1)),
                pc=pc_range
            )
        output:
            [
                "results/bardensr/decoded/decoded_" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_l1"
                + str(l1) + "_peak" + str(peak) + ".csv"
                for l1, peak in bardensr_combos
            ]
        params:
            sigma=config["sigma"],
            iterations=config["bardensr_iterations"],
            poolsize=config["bardensr_poolsize"],
            combos=bardensr_combos,
            code=code
        wildcard_constraints:
            rep=r"\d+",
            code="|".join(codes_to_simulate["code"])
        #benchmark: "results/benchmarks/bardensr/" + code + "_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
        conda: "../envs/bardensr.yaml"
        script: "../scripts/run_bardensr.py"
