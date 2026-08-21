import itertools

bardensr_combos = list(itertools.product(config["bardensr_l1_penalty"], config["bardensr_peak_thresh"]))
istdeco_combos = list(itertools.product(config["istdeco_threshold_Q"], config["istdeco_threshold_percentile"]))
spotsfirst_combos = list(itertools.product(config["lf"], config["drop_correction"]))


rule score_bardensr:
    # Scores every (l1, peak) decoded file for this condition and picks the best in one job,
    # instead of a separate score_bardensr_params job per combo (each of which just did a cheap
    # pandas groupby) followed by a separate aggregator job.
    input:
        decoded=[
            "results/bardensr/decoded/decoded_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_l1"
            + str(l1) + "_peak" + str(peak) + ".csv"
            for l1, peak in bardensr_combos
        ],
        truth="results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
        codebook="resources/codebooks/{code}cb.csv"
    output:
        "results/scores/bardensr/score_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        method="bardensr",
        combos=bardensr_combos
    wildcard_constraints:
        rep=r"\d+",
        nbarcodes=r"\d+",
        nnonspec=r"\d+",
        pdrop=r"[\d.]+",
        rstdv=r"[\d.]+"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_bardensr_best.py"


rule score_istdeco:
    # Scores every (Q, tau) decoded file for this condition and picks the best in one job,
    # instead of a separate score_istdeco_params job per combo followed by a separate
    # aggregator job.
    input:
        decoded=[
            "results/istdeco/decoded/decoded_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_Q"
            + str(Q) + "_tau" + str(tau) + ".csv"
            for Q, tau in istdeco_combos
        ],
        truth="results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
        codebook="resources/codebooks/{code}cb.csv"
    output:
        "results/scores/istdeco/score_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        method="istdeco",
        combos=istdeco_combos
    wildcard_constraints:
        rep=r"\d+",
        nbarcodes=r"\d+",
        nnonspec=r"\d+",
        pdrop=r"[\d.]+",
        rstdv=r"[\d.]+"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_istdeco_best.py"


rule score_graph_iss_params:
    input:
        decoded="results/graph_iss/decoded/decoded_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_k1{k1}.csv",
        truth="results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
        codebook="resources/codebooks/{code}cb.csv"
    output:
        "results/scores/graph_iss/score_params_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_k1{k1}.csv"
    params:
        method="graph_iss"
    wildcard_constraints:
        rep=r"\d+",
        nbarcodes=r"\d+",
        nnonspec=r"\d+",
        pdrop=r"[\d.]+",
        rstdv=r"[\d.]+",
        k1=r"[\d.]+"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_decoded.py"


rule score_graph_iss:
    input:
        expand(
            "results/scores/graph_iss/score_params_{{code}}_rep_{{rep}}_nbarcodes{{nbarcodes}}_nnonspec{{nnonspec}}_pdrop{{pdrop}}_rstdv{{rstdv}}_k1{k1}.csv",
            k1=config["graph_iss_k1"]
        )
    output:
        "results/scores/graph_iss/score_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        method="graph_iss"
    wildcard_constraints:
        rep=r"\d+",
        nbarcodes=r"\d+",
        nnonspec=r"\d+",
        pdrop=r"[\d.]+",
        rstdv=r"[\d.]+"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_graph_iss_best.py"


rule score_spotsfirst:
    # Adds truth-based scoring to every (lf, drc) sum_stats/decoded pair for this condition and
    # picks the best in one job, instead of a separate spotsfirst_add_truth job per combo
    # followed by a separate aggregator job.
    input:
        sum_stats=[
            "results/spotsfirst/sum_stats/sum_stats_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_lf"
            + str(lf) + "_drc" + str(drc) + ".csv"
            for lf, drc in spotsfirst_combos
        ],
        decoded=[
            "results/spotsfirst/decoded/decoded_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_lf"
            + str(lf) + "_drc" + str(drc) + ".csv"
            for lf, drc in spotsfirst_combos
        ],
        truth="results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
        codebook="resources/codebooks/{code}cb.csv"
    output:
        "results/scores/spotsfirst/score_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        method="Untangle: spots first"
    wildcard_constraints:
        rep=r"\d+",
        nbarcodes=r"\d+",
        nnonspec=r"\d+",
        pdrop=r"[\d.]+",
        rstdv=r"[\d.]+"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_spotsfirst_best.py"


rule score_compressed_sensing:
    input:
        decoded_perlambda="results/compressed_sensing/decoded_perlambda/nc_perlambda_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        truth="results/truth/ground_truth_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_sim.csv",
        codebook="resources/codebooks/{code}cb.csv"
    output:
        "results/scores/compressed_sensing/score_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    params:
        method="Untangle: CS l0"
    conda: "../envs/python_base.yaml"
    script: "../scripts/score_cs_bestlambda.py"
