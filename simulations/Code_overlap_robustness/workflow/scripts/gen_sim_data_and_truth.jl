#using Images
#using FileIO
using DataFrames
using CSV
using DelimitedFiles
using Distributions

"""
    sim_true(n)

Simulate n objects in a 1 x 1 field of view. The species of each object is a
random integer from 1 to 1267. It's the location is drawn from the uniform
distribution in the field of view. Returns a data frame of each object's species
x, and y position.
"""

sim_true(n, max_species) = DataFrame(gene = rand(1:max_species, n), x = rand(n), y = rand(n))

function get_n_q_w(cb)
    ncws, n = size(cb)
    q = length(unique(cb))
    if typeof(cb[1,1]) == String7 || typeof(cb[1,1]) == String
        w = maximum(sum(String.(cb) .!= "0", dims=2))
    else
        w = maximum(sum(.~ iszero.(cb), dims=2))
    end
    [n, q, w]
end

savestring_2_pseudocolor_df = DataFrame(CSV.File(snakemake.input[2]))
savestring_2_pseudocolor = Dict((String(row.ffelem), UInt8(eval(row.pseudocolor))) for row in eachrow(savestring_2_pseudocolor_df))
"""
    encode(true_locs :: DataFrame, cb :: Array{Int8,2})

Given a data frame of true locations, and a codebook, return a dataframe giving
the hybridization of the encoded dots, their gene_number, their true locations, and
coefficient.
"""
function encode(true_locs :: DataFrame, cb, pdrop)
    n, q, wt = get_n_q_w(cb)
    n = UInt8(n)
    q = UInt8(q)
    wt = UInt8(wt)
    println("n: $n, q: $q, w: $w")
    #n = UInt8(length(cb[1,:]))
    #wt = sum(cb[1,:] .!= 0)
    #q = UInt8(length(unique(cb)))

    qm1 = UInt8(q-1)

    npts = nrow(true_locs) #length(true_locs.species)
    dot = Array(1:(npts*wt))
    species = Int64[]
    x = Float64[]
    y = Float64[]
    #hyb = Int8[]
    hyb = UInt8[]
    blocks = UInt8[]
    pseudocolors = UInt8[]
    
    for i = 1:npts
        target = true_locs[i, :]
        codeword = cb[target.gene,:]
        for block = 0x01:n
            if codeword[block] != 0 && n > wt && codeword[block] != "0"
                push!(species, target.gene)
                push!(x, target.x)
                push!(y, target.y)
                #dot_hyb = qm1*(block-1)+codeword[block]
                push!(blocks, block)
                if q == 8 || q == 9 
                    pcolor = savestring_2_pseudocolor[String(codeword[block])]
                    push!(pseudocolors, pcolor) #UInt8(eval(savestring_2_pseudocolor[String(codeword[block])])))
                else
                    push!(pseudocolors, codeword[block])
                end
                #push!(hyb, UInt8(dot_hyb))
            elseif n == wt
                push!(species, target.gene)
                push!(x, target.x)
                push!(y, target.y)
                push!(blocks, block)
                #dot_hyb = q*(block-1)+codeword[block]
                #codeword[pos] == 0x00 ? dot_hyb += 0x10 : nothing
                if codeword[block] == 0
                    push!(pseudocolors, UInt8(q))
                else
                    push!(pseudocolors, codeword[block])
                end
                #codeword[block] == 0x00 ? dot_hyb += UInt8(q) : nothing
                #push!(hyb, UInt8(dot_hyb))
            end
        end
    end
    #if n > wt
    pnts = DataFrame(block=blocks,pseudocolor=pseudocolors,cellid=1,gene=species,x=x,y=y)
    #else
    #    pnts = DataFrame(hyb=hyb,gene=species,x=x,y=y)
    #end

    # drop dots
    to_drop = collect(1:nrow(pnts)) # rand(1:wt , npts) + wt*Array(0:(npts-1))
    to_drop = to_drop[rand(nrow(pnts)) .< pdrop]
    deleteat!(pnts, to_drop)

    pnts[!, "z"] .= 1.0
    pnts[!, "w"] .= 1.0
    pnts[!, "s"] .= 1.0
    pnts
end


"""
    draw_rand_dots(n)

Draw random off target dots that simulate probes sticking to some unintended
object that is not of interest.
"""
draw_rand_dots(n_rand_dots, n_block, q) = DataFrame(gene=0, cellid=1, block=rand(1:n_block,n_rand_dots), pseudocolor=rand(1:(q-1),n_rand_dots), x=rand(n_rand_dots), y=rand(n_rand_dots), z=ones(n_rand_dots), w =ones(n_rand_dots),s=ones(n_rand_dots))

#add_rand_dots(pnts, n_rand, n_block, q) = sort(vcat(pnts, draw_rand_dots(n_rand, n_block, q)),[:block, :pseudocolor])

function add_rand_dots(pnts, p, n_block, q, zeros_probed)
    if zeros_probed
        start_pc =  0
        npc = q 
    else
        start_pc = 1
        npc = q-1
    end

    rand_dots = DataFrame(
                            gene=0,
                            cellid=1, 
                            block=vcat(repeat.(collect.(1:n_block), npc)...), 
                            pseudocolor=repeat(start_pc:(q-1), n_block),
                            x=rand(n_block*npc), 
                            y=rand(n_block*npc), 
                            z=ones(n_block*npc), 
                            w =ones(n_block*npc),
                            s=ones(n_block*npc)
    )
    rand_dots = rand_dots[rand(Float64, nrow(rand_dots)) .<= p, :]
    println("nrand_dots ", nrow(rand_dots))
    return sort(vcat(pnts, rand_dots))
end


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
#cb = readdlm(snakemake.input[1], ',', UInt8)

cb = DataFrame(CSV.File(cb))
println("read codebook")
if typeof(cb.block1[1]) == String7 || typeof(cb.block1[1]) == String
    cb = Matrix((String.(cb[!, 2:end])))
else
    cb = Matrix(UInt8.(cb[!, 2:end]))
end
#cb = cb .% maximum(cb)
n_symbols,q,w = get_n_q_w(cb)
n_cws, n_symbols = size(cb)

fov_width = snakemake.params["fov_width"]
rstdv =  snakemake.params["rstdv"]
sigma = snakemake.params["sigma"]
ntargets = parse(Int64, snakemake.wildcards["nbarcodes"])
prand_dot  = parse(Float64, snakemake.wildcards["prand"])
println("prand_dot: $prand_dot")

code = snakemake.wildcards["code"]
zeros_probed = code[1:7] == "seqFISH"

pdrop = parse(Float64, snakemake.wildcards["pdrop"])

true_locs = sim_true(ntargets, n_cws)#, fov_width)

CSV.write(snakemake.output[1], true_locs)
println("encoding")
pnts = encode(true_locs, cb, pdrop)
pnts[!,"s"] .= sigma 
#println(pnts)
add_localization_errors!(pnts, rstdv)
pnts = add_rand_dots(pnts, prand_dot, n_symbols, q, zeros_probed)


CSV.write(snakemake.output[2], pnts)
println("done")
"""
gblur = GaussBlur2D(sigma, sigma, fov_width)

for hyb in 1:80
    fname = snakemake.output[hyb+2]
    hyb_pnts = Array(Array(pnts[pnts.hyb .== hyb, ["x", "y", "s", "w"]])')
    hyb_im = reshape(phi(gblur, hyb_pnts), fov_width, fov_width)
    save(fname, block.(UInt16, hyb_im))
end


"""