rule spotsfirst_find_barcode_candidates:
    input:
        "results/dots/dots_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        "resources/codebooks/{code}cb.csv",
        "resources/codebooks/{code}_H.csv",
        "workflow/envs/cs_julia_environment/",
        "results/cs_julia_environment_installed.txt"
    params:
        lf=min(config["lf"]),
        drc=max(config["drop_correction"]),
        lat_thresh=config["lat_thresh"]
    output:
        "results/spotsfirst/candidate_barcodes/candidate_barcodes_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv"
    benchmark: "results/benchmarks/spotsfirst_candidates/{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.tsv"
    # Shares a group id with spotsfirst_choose_barcodes and score_spotsfirst for this same
    # condition so cluster execution submits them as one scheduler job instead of several.
    group: "spotsfirst_{code}_{rep}_{nbarcodes}_{nnonspec}_{pdrop}_{rstdv}"
    script: "../scripts/find_save_barcode_candidates.jl"


rule spotsfirst_choose_barcodes:
    input:
        "results/spotsfirst/candidate_barcodes/candidate_barcodes_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        "results/dots/dots_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}.csv",
        "resources/codebooks/{code}cb.csv",
        "resources/codebooks/{code}_H.csv",
        "workflow/envs/cs_julia_environment/",
        "results/cs_julia_environment_installed.txt"
    params:
        skip_thresh="auto",
        skip_density_thresh=100,
        lat_thresh=config["lat_thresh"]
    output:
        "results/spotsfirst/decoded/decoded_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_lf{lf}_drc{drc}.csv",
        "results/spotsfirst/sum_stats/sum_stats_{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_lf{lf}_drc{drc}.csv",
    benchmark: "results/benchmarks/spotsfirst_ilp/{code}_rep_{rep}_nbarcodes{nbarcodes}_nnonspec{nnonspec}_pdrop{pdrop}_rstdv{rstdv}_lf{lf}_drc{drc}.tsv"
    # Same group id as spotsfirst_find_barcode_candidates for this condition (deliberately
    # omits lf/drc so every (lf, drc) instance joins the one group job for its condition).
    group: "spotsfirst_{code}_{rep}_{nbarcodes}_{nnonspec}_{pdrop}_{rstdv}"
    script: "../scripts/choose_barcodes_from_saved_candidates.jl"
