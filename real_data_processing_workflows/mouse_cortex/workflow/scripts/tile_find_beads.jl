using Pkg
Pkg.activate(snakemake.input[2])

using SeqFISH_ADCG
using DelimitedFiles
using CSV
using Statistics
using Images, FileIO

#parameters
sigma_lb = snakemake.params["sigma_lb"]
sigma_ub = snakemake.params["sigma_ub"]
tau = 2.0*10^12
final_loss_improvement = snakemake.params["final_loss_improvement"]
min_weight = snakemake.params["min_weight"]
max_iters = snakemake.params["max_iters"]
max_cd_iters = snakemake.params["max_cd_iters"]
min_allowed_separation = Float64(snakemake.params["min_allowed_separation"])


img = load(snakemake.input[1])
img = reinterpret.(UInt16, channelview(img))

threshold = exp(mean(log.(img)))

points_w_dup, records = fit_2048x2048_img_tiles(img, sigma_lb, sigma_ub, tau, final_loss_improvement, min_weight, max_iters, max_cd_iters, threshold)#mean_int)

points_1 = remove_duplicates(points_w_dup, img, sigma_lb, sigma_ub, tau, threshold, min_allowed_separation)
points_2 = remove_duplicates(points_1, img, sigma_lb, sigma_ub, tau, threshold, min_allowed_separation)
points = remove_duplicates(points_2, img, sigma_lb, sigma_ub, tau, threshold, min_allowed_separation)



CSV.write(snakemake.output[1], points)

