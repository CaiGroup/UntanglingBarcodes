using Pkg
Pkg.activate(snakemake.input[3])

using SeqFISH_ADCG
using Images
using FileIO
using DataFrames
using CSV
using DelimitedFiles
using Distributions


function get_n_q_w(cb :: Matrix)
    ncws, n = size(cb)
    q = length(unique(cb))
    if typeof(cb[1,1]) <: AbstractString
        w = maximum(sum(String.(cb) .!= "0", dims=2))
    else
        w = maximum(sum(.~ iszero.(cb), dims=2))
    end
    [n, q, w]
end

"""
    encode(true_locs :: DataFrame, cb :: Array{Int8,2})

Given a data frame of true locations, and a codebook, return a dataframe giving
the hybridization of the encoded dots, their gene_number, their true locations, and
coefficient.
"""
function encode(true_locs :: DataFrame, cb :: DataFrame, pbind_secondary, pdrop)
    n, q, wt = get_n_q_w(Matrix(select(cb, Not(:gene))))
    n = UInt8(n)
    q = UInt8(q)
    wt = UInt8(wt)

    qm1 = UInt8(q-1)

    npts = nrow(true_locs)
    dot = Array(1:(npts*wt))
    species = Int[]
    primary = Int[]
    x = Float64[]
    y = Float64[]
    #hyb = Int8[]
    hyb = UInt8[]
    blocks = UInt8[]
    pseudocolors = UInt8[]

    for i = 1:npts
        target = true_locs[i, :]
        codeword = Matrix(select(cb, Not(:gene))[target["gene"] .== cb.gene,:])
        codeword = reshape(codeword, n)
        for block = 0x01:n
            if codeword[block] != 0 && n > wt && codeword[block] != "0"
                push!(species, target.gene)
                push!(x, target.x)
                push!(y, target.y)
                push!(blocks, block)
                push!(primary, target.primary)
                push!(pseudocolors, codeword[block])
            elseif n == wt
                push!(species, target.gene)
                push!(x, target.x)
                push!(y, target.y)
                push!(primary, target.primary)
                push!(pseudocolors, codeword[block])
                push!(blocks, block)
                #codeword[block] == 0x00 ? dot_hyb += UInt8(q) : nothing
            end
        end
    end
    #if n > wt
    pnts = DataFrame(block=blocks,pseudocolor=pseudocolors,cellid=1,gene=species,x=x,y=y,primary=primary)
    #else
    #    pnts = DataFrame(hyb=hyb,gene=species,x=x,y=y,primary=primary)
    #end

    pnts[!, "z"] .= 1.0
    nsecondary = [rand(Distributions.Binomial(nprim, pbind_secondary)) for nprim in pnts.primary]
    pnts[!, "w"] .= snakemake.params["dot_intensity"] *nsecondary/20#10000
    pnts[!, "s"] .= 1.0
    pnts
end


"""
    draw_rand_dots(n)

Draw random off target dots that simulate probes sticking to some unintended
object that is not of interest.
"""
pbind_prim = snakemake.params["pbind_primary"]
draw_rand_dots(n_rand_dots, n_block, q) = DataFrame(gene=0, cellid=1, block=rand(1:n_block, n_rand_dots), pseudocolor=rand(1:(q-1), n_rand_dots), x=rand(n_rand_dots), y=rand(n_rand_dots), z=ones(n_rand_dots), w =ones(n_rand_dots),s=ones(n_rand_dots),primary=rand(Distributions.Binomial(20,pbind_prim), n_rand_dots))

add_rand_dots(pnts, n_rand, n_block, q) = sort(vcat(pnts, draw_rand_dots(n_rand, n_block, q)),[:block, :pseudocolor])


"""
    draw_localization_errors(n, rstdv)

Draw vectors of x and y localization errors for n points.
"""
function draw_localization_errors(n, rstdv)
    r_error = rand(Normal(0,rstdv), n)
    θ = 2*pi*rand(Float64, n)
    x_error = r_error.*sin.(θ)
    y_error = r_error.*cos.(θ)
    return x_error, y_error
end

"""
    add_localization_errors(pnts :: DataFrame, x_errors, y_errors)
"""
function add_localization_errors!(pnts :: DataFrame, x_errors, y_errors)
    pnts[:,"x"] .= pnts[:,"x"] + x_errors
    pnts[:,"y"] .= pnts[:,"y"] + y_errors
end

"""
    add_localization_errors(pnts :: DataFrame)
"""
function add_localization_errors!(pnts :: DataFrame, rstdv)
    x_errors, y_errors = draw_localization_errors(length(pnts.x), rstdv)
    add_localization_errors!(pnts, x_errors, y_errors)
end

cb = snakemake.input[1]

cb = DataFrame(CSV.File(cb))
println("read codebook")
n_symbols,q,w = get_n_q_w(Matrix(select(cb, Not(:gene))))

fov_width = snakemake.params["fov_width"]
rstdv =  snakemake.params["rstdv"]
sigma = snakemake.params["sigma"]
n_rand_dots = parse(Int64, snakemake.wildcards["nnonspec"])
pdrop = parse(Float64, snakemake.wildcards["pdrop"])
pbind_secondary = snakemake.params["pbind_secondary"]

true_locs = DataFrame(CSV.File(snakemake.input[2]))

CSV.write(snakemake.output[1], true_locs)
println("encoding")
pnts = encode(true_locs, cb, pbind_secondary, pdrop)
pnts[!,"s"] .= sigma
add_localization_errors!(pnts, rstdv)
pnts = add_rand_dots(pnts, n_rand_dots, n_symbols, q)


CSV.write(snakemake.output[2], pnts)
println("done")

gblur = GaussBlur2D(sigma, sigma, fov_width)

for i in 3:length(snakemake.output)
    fname = snakemake.output[i]
    println("saving $fname")
    r = parse(Int, split(fname, "_")[4])

    pseudocolor = parse(Int, split(fname, "_")[6])

    hyb_pnts = Array(Array(pnts[(pnts.block .== r) .&& (pnts.pseudocolor .== pseudocolor), ["x", "y", "s", "w"]])')
    hyb_im = reshape(phi(gblur, hyb_pnts), fov_width, fov_width)
    save(fname, round.(UInt16, hyb_im))
end