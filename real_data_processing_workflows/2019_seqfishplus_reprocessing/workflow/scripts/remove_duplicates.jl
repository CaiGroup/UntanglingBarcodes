using DelimitedFiles
using CSV
using DataFrames
using Statistics
using SeqFISH_ADCG #SparseInverseProblems
using Images, FileIO


ch = snakemake.wildcards["ch"]

sigma_ub = Float64(snakemake.params["sigma_ub"][ch])
sigma_lb = Float64(snakemake.params["sigma_lb"][ch])
tau = 2.0*10^12
min_allowed_separation = Float64(snakemake.params["min_allowed_separation"])

img = load(snakemake.input[1])
img = reinterpret.(UInt16, channelview(img))
#img = readdlm(snakemake.input[1])

points_w_dup = DataFrame(CSV.File(snakemake.input[2]))

labeled_img = load(snakemake.input[3])
if typeof(channelview(labeled_img)[1,1]) == Normed{UInt16,16}
    labeled_img = reinterpret.(UInt16, channelview(labeled_img))
elseif typeof(channelview(labeled_img)[1,1]) == Normed{UInt8,8}
    labeled_img = reinterpret.(UInt8, channelview(labeled_img)) #if their are fewer than 255 cells...
else
    error("Pixels in labeled image are neigher UInt16s or UInt8s")
end
noise_mean = 0.0

points_rm1 = remove_duplicates(points_w_dup, img, sigma_lb, sigma_ub, tau, noise_mean, min_allowed_separation)
points_rm2 = remove_duplicates(points_rm1, img, sigma_lb, sigma_ub, tau, noise_mean, min_allowed_separation)
points = remove_duplicates(points_rm2, img, sigma_lb, sigma_ub, tau, noise_mean, min_allowed_separation)



xr = round.(Int64, points.x)
yr = round.(Int64, points.y)

if size(points)[1] == 0
    pos = []
    ch = []
    hyb = []
else
    ysize, xsize = size(labeled_img)
    xr[xr .< 1] .= 1
    xr[xr .> xsize] .= xsize
    yr[yr .< 1] .= 1
    yr[yr .> ysize] .= ysize

    pos = snakemake.wildcards["pos"]
    ch = snakemake.wildcards["ch"]
    hyb = snakemake.wildcards["hyb"]
end


cellid = [labeled_img[yr[i], xr[i]] for i in 1:length(xr)]
points[!,"pos"] .= pos
points[!,"ch"] .= ch
points[!,"hyb"] .= hyb
points[!,"cellid"] = cellid



CSV.write(snakemake.output[1], points)
