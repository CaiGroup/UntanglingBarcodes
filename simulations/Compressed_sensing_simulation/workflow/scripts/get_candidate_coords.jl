using Pkg
Pkg.activate(snakemake.input[2])

using Images
using FileIO
using DataFrames
using CSV
using Plots
using DelimitedFiles
using MultivariateStats
using GLMNet

i = 1
fov_width = snakemake.params["fov_width"]
grid_size = snakemake.params["grid_size"]
min_intensity = snakemake.params["min_intensity"]


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

cb = DataFrame(CSV.File(snakemake.input[1]))

n, q, w = get_n_q_w(Matrix(select(cb,Not(:gene))))

hyperstack = zeros(n*(q-1),1,fov_width,fov_width)
rounds = []
pseudocolors = []


for i in 4:length(snakemake.input)
    im_fname = snakemake.input[i]
    r = parse(Int64, split(im_fname, "_")[4])
    pc = parse(Int64, split(im_fname, "_")[6])
    s = (q-1)*(r -1) + pc
    push!(rounds, r)
    push!(pseudocolors, pc)
    img = reinterpret.(UInt16, channelview(load(im_fname)))
    hyperstack[s,1,:,:] .= img
end

psf(x,y,z) = exp(-(x^2 + y^2 + z^2)/2.0)


psfs = fill(psf, n*(q-1))
psf_depths = ones(n*(q-1))
psf_heights = 6 .* ones(n*(q-1))
psf_widths = 6 .* ones(n*(q-1))
min_intensities = 8000.0 .* ones(n*(q-1))

lat_search_radius = 2
z_search_radius = 1
zeros_probed = false


function get_cand_dots(hyperstack, psfs, psf_depths, psf_heights, psf_widths, min_intensities, rounds, pseudocolors)

    conv_threshold = get_conv_threshold.(psfs, psf_depths, psf_heights, psf_widths, min_intensities)

    im_cand_dot_dfs = []
    for i in 1:length(psfs)
        println("i $i")

        #covolve image and psf
        im_candidate_dot_coords = get_cand_dot_coords(hyperstack[i, :,:,:], psfs[i], min_intensities[i], conv_threshold[i], psf_heights[i], psf_widths[i], psf_depths[i])
        
        #add pseudocolor and round columns
        im_candidate_dot_coords[!, "block"] .= rounds[i]
        im_candidate_dot_coords[!, "pseudocolor"] .= pseudocolors[i]
        
        push!(im_cand_dot_dfs, im_candidate_dot_coords)

    end

    println(nrow.(im_cand_dot_dfs))
    candidate_dot_coords = vcat(im_cand_dot_dfs...)
    return candidate_dot_coords
end

get_conv_threshold(psf, depth, height, width, min_intensity) = min_intensity *sum([psf(x -width/2, y - height/2, z - depth/2) for x in 1:width, y in 1:height, z in 1:depth] .^2)


function get_cand_dot_coords(imagestack, psf, min_intensity, thresh, psf_height, psf_width, psf_depth)
    im_depth, im_height, im_width = size(imagestack)
    cand_dot_coords = []
    observed = zeros(im_depth*im_height*im_width)
    psf_mod = zeros(im_depth*im_height*im_width)
    for z1 in 1:0.5:im_depth, y1 in 1:grid_size:im_height, x1 in 1:grid_size:im_width
        conv_val = 0
        observed .= 0.0
        psf_mod .= 0.0
        i = 1
        for z2 in maximum([1, floor(Int64, z1 - psf_depth/2)]):minimum([im_depth, ceil(Int64, z1 + psf_depth/2)])
            for y2 in maximum([1, floor(Int64, y1 - psf_height/2)]):minimum([im_height, ceil(Int64, y1 + psf_height/2)])
                for x2 in maximum([1, floor(Int64, x1 - psf_width/2)]):minimum([im_width, ceil(Int64, x1 + psf_width/2)])
                    conv_val += imagestack[z2, y2, x2] * psf(x1-x2, y2-y2, z1-z2)
                    observed[i] = imagestack[z2, y2, x2]
                    psf_mod[i] = psf(x1-x2, y2-y2, z1-z2)
                end
            end
        end
        if conv_val > thresh

            push!(cand_dot_coords, [z1, y1, x1])
        end
    end
    if length(cand_dot_coords) > 0
        df = DataFrame(hcat(cand_dot_coords...)', [:z, :y, :x])
        df = lasso_cand_dots(imagestack, psf, df, psf_width, psf_height, psf_depth)
        df[!, "s"] .= 0
        df[!, "s_xy"] .= 0
        df[!, "s_z"] .= 0
        df[!, "w"] .= 1
    else
        df = DataFrame()
        df[!,"z"] = []
        df[!,"y"] = []
        df[!,"x"] = []
        df[!,"s"] = []
        df[!,"s_xy"] = []
        df[!,"s_z"] = []
        df[!,"w"] = []
    end
    return df
end


function lasso_cand_dots(img, psf, cand_dots, psf_width, psf_height, psf_depth)
    im_depth, im_height, im_width = size(img)

    npixels = length(img)
    observed = zeros(npixels)
    A = zeros(npixels, nrow(cand_dots))

    i = 1
    for z in 1:im_depth, y in 1:im_height, x in 1:im_width
        observed[i] =  img[z,y,x]
        for j in 1:nrow(cand_dots)
            if (abs(x - cand_dots.x[j]) <= psf_width/2) && (abs(y - cand_dots.y[j]) <= psf_height/2) && (abs(z - cand_dots.z[j]) <= psf_depth/2)
                A[i,j] = psf(x - cand_dots.x[j], y - cand_dots.y[j], z - cand_dots.z[j])
            end
        end
        i += 1
    end

    constraints = zeros(2, nrow(cand_dots))
    constraints[2,:] .= 2^16
    lasso_mod = glmnet(A, observed, constraints = constraints)
    lassoed_candidates = cand_dots[findall(c -> c > 0, reshape(sum(lasso_mod.betas, dims=2), nrow(cand_dots))),:]
    return lassoed_candidates
end

candidate_dot_coords = get_cand_dots(hyperstack, psfs, psf_depths, psf_heights, psf_widths, min_intensities, rounds, pseudocolors)

CSV.write(snakemake.output[1], candidate_dot_coords)