#using NetCDF
#using DelimitedFiles
using CSV
using Images, FileIO
using SeqFISH_ADCG #SparseInverseProblems
using DataFrames

pos = snakemake.wildcards["pos"]
ch = snakemake.wildcards["ch"]
hyb = snakemake.wildcards["hyb"]

#parameters
sigma_lb = snakemake.params["sigma_lb"][ch]
sigma_ub = snakemake.params["sigma_ub"][ch]
tau = 2.0*10^12
final_loss_improvement = snakemake.params["final_loss_improvement"][ch]
min_weight = snakemake.params["min_weight"][ch]
max_iters = snakemake.params["max_iters"]
max_cd_iters = snakemake.params["max_cd_iters"]

#img = ncread(snakemake.input[1], "fluorescence")
# img = load("HybCycle_9_ch_640_pos_5.tif")
img = load(snakemake.input[1])
#println("img: ", typeof(img))
img = reinterpret.(UInt16, channelview(img))
#println("img reinterpreted:", typeof(img))

#labeled_img = load("hyb_9_pos_5.tif")#snakemake.input[2])
labeled_img = load(snakemake.input[2])
if typeof(channelview(labeled_img)[1,1]) == Normed{UInt16,16}
    labeled_img = reinterpret.(UInt16, channelview(labeled_img))
elseif typeof(channelview(labeled_img)[1,1]) == Normed{UInt8,8}
    labeled_img = reinterpret.(UInt8, channelview(labeled_img)) #if their are fewer than 255 cells...
else
    #println("they are type: ", typeof(channelview(labeled_img)[1,1]))
    error("Pixels in labeled image are neither UInt16s or UInt8s")# they are type: ", typeof(channelview(labeled_img)[1,1])")
    #println("they are type: ", typeof(channelview(labeled_img)[1,1]))
end

#mask = labeled_img .> 0
#asked_img = img .* mask

threshold = 0.0

#points = fit_2048x2048_img_tiles(masked_img, sigma_lb, sigma_ub, tau, final_loss_improvement, min_weight, max_iters, max_cd_iters, threshold)
points_w_dup, records = fit_2048x2048_img_tiles(img, sigma_lb, sigma_ub, tau, final_loss_improvement, min_weight, max_iters, max_cd_iters, threshold)


xr = round.(Int64, points_w_dup.x)
yr = round.(Int64, points_w_dup.y)

#ysize, xsize = size(masked_img)
ysize, xsize = size(img)


xr[xr .< 1] .= 1
xr[xr .> xsize] .= xsize
yr[yr .< 1] .= 1
yr[yr .> ysize] .= ysize

cellid = [labeled_img[yr[i], xr[i]] for i in 1:length(xr)]
points_w_dup[!,"pos"] .= pos
points_w_dup[!,"ch"] .= ch
points_w_dup[!,"hyb"] .= hyb
points_w_dup[!,"cellid"] = cellid
points_w_dup = points_w_dup[points_w_dup.cellid .!= 0, :]

CSV.write(snakemake.output[1], points_w_dup)

# https://gitlab.com/jawhitect/sparseinverseproblems
