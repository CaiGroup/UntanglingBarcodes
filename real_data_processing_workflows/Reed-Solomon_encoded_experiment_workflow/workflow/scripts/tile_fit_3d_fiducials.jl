using Images, FileIO, SeqFISH_ADCG, CSV, DataFrames

pos = snakemake.wildcards["pos"]
ch = snakemake.wildcards["ch"]
if "hyb" in keys(snakemake.wildcards)
    hyb = snakemake.wildcards["hyb"]
end

imst_ = load(snakemake.input[1])
imst = reinterpret.(UInt16, channelview(imst_))
#imst = zeros(UInt16, 2048, 2048, 3)
#println("size(imst_): ", size(imst))
#for z in 1:17
#    imst[:,:,z] = imst_[z,:,:]
#end

#ch = 0 #snakemake.wildcards["ch"]
#sigma_lb = snakemake.params["sigma_lb"][ch]
#sigma_ub = snakemake.params["sigma_ub"][ch]

main_tile_width = snakemake.params["tile_main_width"]
tile_overlap = snakemake.params["tile_overlap"]
tile_depth = snakemake.params["tile_depth"]
tile_depth_overhang = snakemake.params["tile_depth_overhang"]
sigma_xy_lb = snakemake.params["sigma_xy_lb"]
sigma_xy_ub = snakemake.params["sigma_xy_ub"]
sigma_z_lb = snakemake.params["sigma_z_lb"]
sigma_z_ub = snakemake.params["sigma_z_ub"]
tau = 2.0*10^12
final_loss_improvement = snakemake.params["final_loss_improvement"]
min_weight = snakemake.params["min_weight"][ch]
max_iters = snakemake.params["max_iters"]
max_cd_iters = snakemake.params["max_cd_iters"]
min_allowed_separation = snakemake.params["min_allowed_separation"]

points_w_dup, records =  fit_img_tiles(imst, main_tile_width, tile_overlap, sigma_xy_lb, sigma_xy_ub,
    tau, final_loss_improvement,min_weight, max_iters, max_cd_iters, 0.0, fit_alg="ADCG")

#points_w_dup, records =  fit_stack_tiles(imst, main_tile_width, tile_overlap, tile_depth, tile_depth_overhang, sigma_xy_lb, sigma_xy_ub,
#    sigma_z_lb, sigma_z_ub, final_loss_improvement,min_weight, max_iters, max_cd_iters, "ADCG")
#points_w_dup =  fit_stack_tiles(imst, main_tile_width, tile_overlap, sigma_xy_lb, sigma_xy_ub,
#    sigma_z_lb, sigma_z_ub, final_loss_improvement,min_weight, max_iters, max_cd_iters)
#points = remove_duplicates3d(points_w_dup, sigma_xy_lb, sigma_xy_ub, min_allowed_separation)
points = remove_duplicates(points_w_dup, imst, sigma_xy_lb, sigma_xy_ub, tau, 0.0, min_allowed_separation)

xr = round.(Int64, points.x)
yr = round.(Int64, points.y)
zr = zeros(nrow(points)) #round.(Int64, points.z)


#ysize, xsize = size(masked_img)
#zsize, ysize, xsize = size(imst)
#ysize, xsize = size(imst)

#xr[xr .< 1] .= 1
#xr[xr .> xsize] .= xsize
#yr[yr .< 1] .= 1
#yr[yr .> ysize] .= ysize
#zr[zr .< 1] .= 1
#zr[zr .> zsize] .= zsize

points[!,"z"] .= zr
points[!,"pos"] .= pos
points[!,"ch"] .= ch


if "hyb" in keys(snakemake.wildcards)
    points[!,"hyb"] .= hyb
end

CSV.write(snakemake.output[1], points)
#CSV.write(snakemake.output[2], records)
