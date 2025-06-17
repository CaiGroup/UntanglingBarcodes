using Images, FileIO, SeqFISH_ADCG, CSV, DataFrames, Statistics

pos = snakemake.wildcards["pos"]
#if "ch" in keys(snakemake.wildcards)
#    global ch = snakemake.wildcards["ch"]
#else
#    global ch = 594
#end
ch = snakemake.wildcards["ch"]
hyb = snakemake.wildcards["hyb"]

imst_ = load(snakemake.input[1])
imst = reinterpret.(UInt16, channelview(imst_))
#imst = zeros(UInt16, 2048, 2048, 3)
#for z in 1:17
#    imst[:,:,z] = imst_[z,:,:]
#end
#ch = 0 #snakemake.wildcards["ch"]
#sigma_lb = snakemake.params["sigma_lb"][ch]
#sigma_ub = snakemake.params["sigma_ub"][ch]
nzeros = filter(pv -> pv != 0, imst)

#min_weight = quantile(nzeros, 0.98) #2)

main_tile_width = snakemake.params["tile_main_width"]
tile_overlap = snakemake.params["tile_overlap"]
tile_depth = snakemake.params["tile_depth"]
tile_depth_overhang = snakemake.params["tile_depth_overhang"]
sigma_xy_lb = snakemake.params["sigma_xy_lb"][ch]
sigma_xy_ub = snakemake.params["sigma_xy_ub"][ch]
sigma_z_lb = snakemake.params["sigma_z_lb"][ch]
sigma_z_ub = snakemake.params["sigma_z_ub"][ch]
tau = 2.0*10^12
final_loss_improvement = snakemake.params["final_loss_improvement"][ch]
min_weight = snakemake.params["min_weight"][ch]
max_iters = snakemake.params["max_iters"][ch]
max_cd_iters = snakemake.params["max_cd_iters"]
min_allowed_separation = snakemake.params["min_allowed_separation"]


"""
labeled_img = load(snakemake.input[2])
if typeof(channelview(labeled_img)[1,1,1]) == Normed{UInt16,16}
    labeled_img = reinterpret.(UInt16, channelview(labeled_img))
elseif typeof(channelview(labeled_img)[1,1,1]) == Normed{UInt8,8}
    labeled_img = reinterpret.(UInt8, channelview(labeled_img)) #if their are fewer than 255 cells...
else
    #println("they are type: ", typeof(channelview(labeled_img)[1,1]))
    error("Pixels in labeled image are neither UInt16s or UInt8s")# they are type: ", typeof(channelview(labeled_img)[1,1])")
    #println("they are type: ", typeof(channelview(labeled_img)[1,1]))
end
"""
points_w_dup, records =  fit_img_tiles(imst, main_tile_width, tile_overlap, sigma_xy_lb, sigma_xy_ub,
    tau, final_loss_improvement,min_weight, max_iters, max_cd_iters, 0.0, fit_alg="ADCG")


#points_w_dup, records =  fit_stack_tiles(imst, main_tile_width, tile_overlap, tile_depth, tile_depth_overhang, sigma_xy_lb, sigma_xy_ub,
#    sigma_z_lb, sigma_z_ub, final_loss_improvement,min_weight, max_iters, max_cd_iters, "ADCG")
#points_w_dup =  fit_stack_tiles(imst, main_tile_width, tile_overlap, sigma_xy_lb, sigma_xy_ub,
#    sigma_z_lb, sigma_z_ub, final_loss_improvement,min_weight, max_iters, max_cd_iters)
points = remove_duplicates(points_w_dup,imst, sigma_xy_lb, sigma_xy_ub, tau, 0.0, min_allowed_separation)
#points = remove_duplicates3d(points_w_dup, sigma_xy_lb, sigma_xy_ub, min_allowed_separation)

"""
xr = round.(Int64, points.x)
yr = round.(Int64, points.y)
zr = round.(Int64, points.z)


#ysize, xsize = size(masked_img)
zsize, ysize, xsize = size(labeled_img)


xr[xr .< 1] .= 1
xr[xr .> xsize] .= xsize
yr[yr .< 1] .= 1
yr[yr .> ysize] .= ysize
zr[zr .< 1] .= 1
zr[zr .> zsize] .= zsize

cellid = [labeled_img[zr[i], yr[i], xr[i]] for i in 1:length(xr)]
"""     

points[!,"z"] .= zeros(nrow(points))

points[!,"pos"] .= pos
points[!,"ch"] .= ch
points[!,"hyb"] .= hyb
#points[!,"cellid"] = cellid
#points = points[points.cellid .!= 0, :]
println("saving points. nrow: ", nrow(points))
CSV.write(snakemake.output[1], points)
#CSV.write(snakemake.output[2], records)
