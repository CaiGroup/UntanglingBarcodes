using SeqFISH_ADCG
using Images
using FileIO
using DataFrames
using CSV
using DelimitedFiles
using Distributions


"""
    sim_true(n)

Simulate n objects in a fov_width x fov_width field of view. The gene_number of each object is a
random integer from 1 to max_gene_number. It's the location is drawn from the uniform
distribution in the field of view. Returns a data frame of each object's gene_number
x, and y position.
"""
sim_true(n :: Integer, max_gene_number ::Integer, fov_width :: Integer) = DataFrame(gene_number = rand(1:max_gene_number, n), x = rand(n) .* fov_width, y = rand(n) .* fov_width)

function get_n_q_w(cb)
    ncws, n = size(cb)
    q = length(unique(cb))
    cw_nonzeros = [sum(cb[i, :] .!= -) for i in ncws]
    w =  maximum(cw_nonzeros)
    #w = sum(cb[1,:] .!= 0)
    [n, q, w]
end

"""
    encode(true_locs :: DataFrame, cb :: Array{Int8,2})

Given a data frame of true locations, and a codebook, return a dataframe giving
the hybridization of the encoded dots, their gene_number, their true locations, and
coefficient.
"""
function encode(true_locs :: DataFrame, cb)
    n, q, wt = get_n_q_w(cb)
    n = UInt8(n)
    q = UInt8(q)
    wt = UInt8(wt)
    #n = UInt8(length(cb[1,:]))
    #wt = sum(cb[1,:] .!= 0)
    #q = UInt8(length(unique(cb)))

    qm1 = UInt8(q-1)

    npts = nrow(true_locs) #length(true_locs.gene_number)
    #dot = Array(1:(npts*n))
    gene_number = Int64[]
    x = Float64[]
    y = Float64[]
    #hyb = Int8[]
    hyb = UInt8[]
    barcodeid = []
    for i = 1:npts
        target = true_locs[i, :]
        codeword = cb[target.gene_number,:]
        #for pos = 1:n
        for pos = 0x01:n
            push!(barcodeid, i)
            if codeword[pos] != 0 && n > wt
                push!(gene_number, target.gene_number)
                push!(x, target.x)
                push!(y, target.y)
                dot_hyb = qm1*(pos-1)+codeword[pos]
                push!(hyb, dot_hyb)
                #push!(hyb, UInt8(dot_hyb))
            elseif n == wt
                push!(gene_number, target.gene_number)
                push!(x, target.x)
                push!(y, target.y)
                dot_hyb = q*(pos-1)+codeword[pos]
                #codeword[pos] == 0x00 ? dot_hyb += 0x10 : nothing
                codeword[pos] == 0x00 ? dot_hyb += UInt8(q) : nothing
                push!(hyb, dot_hyb)
                #push!(hyb, UInt8(dot_hyb))
            end
        end
    end
    #pnts = DataFrame(dot_ID=dot,hyb=hyb,gene_number=gene_number,x=x,y=y)
    pnts = DataFrame(hyb=hyb,barcodeid=barcodeid,cellid=1,gene_number=gene_number,x=x,y=y)
    pnts[!, "z"] .= 1.0
    pnts[!, "w"] .= 1000.0
    pnts[!, "s"] .= 1.0
    return pnts
end


"""
    draw_rand_dots(n)

Draw random off target dots that simulate probes sticking to some unintended
object that is not of interest.
"""
draw_rand_dots(n, fov_width) = DataFrame(gene_number=-1, barcodeid =-1, cellid=1, hyb=rand(1:80, n), x=rand(n) .* fov_width, y=rand(n) .* fov_width, z=1.0, w=1.0, s=1.0)

add_rand_dots(pnts, n, fov_width) = sort(vcat(pnts, draw_rand_dots(n, fov_width)), :hyb)

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
cb = readdlm(snakemake.input[1], ',', UInt8)
n_symbols,q,w = get_n_q_w(cb)
n_cws, n_symbols = size(cb)

fov_width = snakemake.params["fov_width"]
rstdv =  snakemake.params["rstdv"]
n_rand_dots = parse(Int64, snakemake.wildcards["nbarcodes"])
sigma = snakemake.params["sigma"]
ntargets = parse(Int64, snakemake.wildcards["nbarcodes"])
n_rand_dots = parse(Int64, snakemake.wildcards["nnonspec"])


true_locs = sim_true(ntargets, n_cws, fov_width)

CSV.write(snakemake.output[1], true_locs)

pnts = encode(true_locs, cb)
pnts[!,"s"] .= sigma 
add_localization_errors!(pnts, rstdv)
pnts = add_rand_dots(pnts, n_rand_dots, fov_width)

CSV.write(snakemake.output[2], pnts)
gblur = GaussBlur2D(sigma, sigma, fov_width)

for hyb in 1:80
    fname = snakemake.output[hyb+2]
    hyb_pnts = Array(Array(pnts[pnts.hyb .== hyb, ["x", "y", "s", "w"]])')
    hyb_im = reshape(phi(gblur, hyb_pnts), fov_width, fov_width)
    save(fname, round.(UInt16, hyb_im))
end


