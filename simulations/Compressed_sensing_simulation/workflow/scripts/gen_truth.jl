using Pkg
Pkg.activate(snakemake.input[2])

using DataFrames
using CSV
using Distributions

function sim_true(n, cb, roi_width, roi_pad, p) 
    sim_genes = rand(1:nrow(cb), n) 
    d = Distributions.Binomial(20,p)

    return DataFrame(gene = sim_genes, x = roi_pad .+ rand(n) .* roi_width, y = roi_pad .+ rand(n) .* roi_width, primary=rand(d, n))
end

cb = DataFrame(CSV.File(snakemake.input[1]))

roi_width = snakemake.params["roi_width"]
roi_pad = snakemake.params["roi_pad"]
pbind_primary = snakemake.params["pbind_primary"]

ntargets = parse(Int64, snakemake.wildcards["nbarcodes"])

true_locs = sim_true(ntargets, cb, roi_width, roi_pad, pbind_primary)

CSV.write(snakemake.output[1], true_locs)
