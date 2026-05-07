using Pkg
using DelimitedFiles

Pkg.activate(snakemake.input[1])
Pkg.instantiate()

open(snakemake.output[1], "w") do io
    writedlm(io, "Done")
end