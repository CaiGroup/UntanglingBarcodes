using SeqFISHPointDecoding
using CSV
using DataFrames
using Statistics
using Profile
using DelimitedFiles

pnts = DataFrame(CSV.File(snakemake.input[1]))#
pnts.hyb = UInt8.(pnts.hyb)
pnts.x = Float64.(pnts.x)
pnts.y = Float64.(pnts.y)
pnts.z = zeros(Float64, nrow(pnts))

H = [1 1 -1 -1;]

#filter!(:hyb=>hyb->hyb != 16 && hyb !=61, pnts)

#lat_thresh = 2.0
z_thresh = 0.0
#code = "E2019"#_neg_ctrl"
code = "E2019_neg_ctrl"
cb_name = snakemake.input[2] #"E2019_cb_all_control.txt"
#cb = readdlm("../codebooks/" * cb_name, UInt8)
cb = readdlm(cb_name, UInt8)

#ndrops = 0#parse(Int64, snakemake.wildcards["dr"])#1
ndrops = parse(Int64, snakemake.wildcards["dr"])#1

pnts.dot_ID = Array(1:nrow(pnts))
ndots = nrow(pnts)

# Cost Parameters
free_dot_cost = 5.0
lat_var_factor = parse(Float64, snakemake.wildcards["lf"])
z_var_factor = 0.0
lw_var_factor = parse(Float64, snakemake.wildcards["wf"])
s_var_factor = parse(Float64, snakemake.wildcards["sf"])
erasure_penalty = 4.0

#lat_thresh = sqrt(20.0/lat_var_factor)*4
#lat_thresh = sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*4 #used for 20210528_results
lat_thresh = sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*3


#Simulated Annealing Parameters
c_final = 1
n_chains = 100
l_chain = 20
cooling_factor = free_dot_cost * 40
mip_sa_thresh = 80
cooling_timescale = exp(cooling_factor/c_final-1)/(n_chains*l_chain)
n_chains = 50
converge_thresh = 100 * ndots


println("cost parameters: $lat_var_factor, $lw_var_factor")


params = SeqFISHPointDecoding.DecodeParams()
set_xy_search_radius(params, lat_thresh)
set_z_search_radius(params, z_thresh)
set_n_allowed_drops(params, ndrops)
set_lat_var_cost_coeff(params, lat_var_factor)
set_z_var_cost_coeff(params, z_var_factor)
set_lw_var_cost_coeff(params, lw_var_factor)
set_s_var_cost_coeff(params, s_var_factor)
set_mip_sa_thresh(params, mip_sa_thresh)
set_n_chains(params, n_chains)
set_l_chains(params, l_chains)
set_converge_thresh(params, converge_thresh)
set_cooling_factor(params, cooling_factor)
set_cooling_timescale(params, cooling_timescale)

mpaths = decode_syndromes!(pnts, cb, H, params)

#@profile mpaths = decode_syndromes!(pnts_cp, lat_thresh, z_thresh, code, ndrops, sa_params)
xs = [pnts.x[dots] for dots in mpaths.cpath]
ys = [pnts.y[dots] for dots in mpaths.cpath]


println("decoded ", sum(pnts.decoded .!= 0), " dots of ", nrow(pnts))
println("Done")

mpaths[!, "xs"] = xs
mpaths[!, "ys"] = ys


println("decoded ", sum(pnts.decoded .!= 0), " dots of ", nrow(pnts))

save_df = pnts[:,append!(Array(1:4), Array(10:11))]

n_mpaths = nrow(mpaths)

n_neg_cntrl_mpaths = sum(mpaths.value .> 3334)
println("Number of mpaths: $n_mpaths")
println("Number of negative control mpaths: $n_neg_cntrl_mpaths")
println("lat var factor: $lat_var_factor")
println("lw var factor: $lw_var_factor")

if cb_name == "Eng2019_647.csv"
    ctrl = "no_neg_ctrl"
else
    ctrl = "w_neg_ctrl"
end

#save_name = "decode_" * ctrl* "_lvf$lat_var_factor"*"_lwvf$lw_var_factor"*"dr$ndrops"*".csv"

sum_stats = DataFrame()
sum_stats["pos"] = parse(Int64, snakemake.wildcards["pos"])
sum_stats["ch"] = parse(Int64, snakemake.wildcards["ch"])
sum_stats["lvf"] = lat_var_factor
sum_stats["lwvf"] = lw_var_factor
sum_stats["svf"] = s_var_factor
sum_stats["dr"] = ndrops
sum_stats["n_barcodes"] = n_mpaths
sum_stats["negative_control_barcodes"] = n_neg_cntrl_mpaths
sum_stats["gene_encoding_barcodes"] = n_mpaths - n_neg_cntrl_mpaths
sum_stats["est_fp_rate"] = (n_neg_cntrl_mpaths*3334/(8000-3334))/n_mpaths


#CSV.write(save_name, save_df)
#CSV.write(snakemake.output[1], save_df)

#CSV.write("mpaths_" * save_name, mpaths)
CSV.write(snakemake.output[1], mpaths)
CSV.write(snakemake.output[2], sum_stats)
CSV.write(snakemake.output[3], save_df)
