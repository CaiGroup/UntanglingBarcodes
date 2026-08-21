using Pkg
Pkg.activate(snakemake.input[5])

using SeqFISHSyndromeDecoding
using CSV
using DataFrames
using Statistics
using DelimitedFiles
using Gurobi

code = snakemake.wildcards["code"]

cpaths_fname = snakemake.input[1]
hyb_pnts_filename =  snakemake.input[2]
cb_name = snakemake.input[3]

pnts = DataFrame(CSV.File(hyb_pnts_filename))
pnts.block = UInt8.(pnts.block)
pnts.pseudocolor = UInt8.(pnts.pseudocolor)
pnts.x = Float64.(pnts.x)
pnts.y = Float64.(pnts.y)
pnts.z = zeros(Float64, nrow(pnts))

H = readdlm(snakemake.input[4], ',')

if any(typeof.(H) .== String7) || any(typeof.(H) .== String) || any(typeof.(H) .== SubString{String})
    H = string.(H)
elseif any(H .< 0)
    H = Int8.(H)
else
    H = UInt8.(H)
end

z_thresh = 0.0
cb = DataFrame(CSV.File(cb_name))

ndrops = parse(Int64, snakemake.wildcards["drc"])

pnts.dot_ID = Array(1:nrow(pnts))
ndots = nrow(pnts)

# Cost Parameters
free_dot_cost = 1.0
lat_var_factor = parse(Float64, snakemake.wildcards["lf"])
z_var_factor = 0.0
lw_var_factor = 0.0 #parse(Float64, snakemake.wildcards["wf"])
s_var_factor = 0.0 #parse(Float64, snakemake.wildcards["sf"])

erasure_penalty = 2*free_dot_cost
lat_thresh = snakemake.params["lat_thresh"] #sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*3


println("cost parameters: $lat_var_factor, $lw_var_factor")

params = DecodeParams()
set_free_dot_cost(params, free_dot_cost)
set_xy_search_radius(params, lat_thresh)
set_z_search_radius(params, z_thresh)
set_n_allowed_drops(params, ndrops)
set_lat_var_cost_coeff(params, lat_var_factor)
set_z_var_cost_coeff(params, z_var_factor)
set_lw_var_cost_coeff(params, lw_var_factor)
set_s_var_cost_coeff(params, s_var_factor)
if code[1:7] == "seqFISH"
    set_zeros_probed(params, true)
else
    set_zeros_probed(params, false)
end
set_skip_thresh(params, 10000000)
set_skip_density_thresh(params, 100000000000)

if snakemake.params["skip_thresh"] == "auto"
    set_skip_thresh(params, nrow(cb))
else
    set_skip_thresh(params, snakemake.params["skip_thresh"])
end
set_skip_density_thresh(params, 100000000000000000000000) #snakemake.params["skip_density_thresh"])

set_erasure_penalty(params, erasure_penalty) #snakemake.params["drop_penalty"])

cpaths_df = DataFrame(CSV.File(cpaths_fname))
cpaths_df[!,"cpath"] = eval.(Meta.parse.(cpaths_df[!,"cpath"]))

if nrow(cpaths_df) == 0
    CSV.write(snakemake.output[1], DataFrame(Dict("cpath"=>[],"cost"=>[],"gene_number"=>[],"x"=>[],"y"=>[],"z"=>[],"cellid"=>[])))
    sum_stats = DataFrame()
    #sum_stats[!, "pos"] = [pos]
    #sum_stats[!, "ch"] = [ch]
    sum_stats[!,"lvf"] = [lat_var_factor]
    sum_stats[!,"lwvf"] = [lw_var_factor]
    sum_stats[!,"svf"] = [s_var_factor]
    sum_stats[!,"dr"] = [ndrops]
    sum_stats[!,"n_barcodes"] = [0]
    sum_stats[!,"negative_control_barcodes"] = [0]
    sum_stats[!,"gene_encoding_barcodes"] = [0]
    sum_stats[!,"est_fp_rate"] = [0]
    CSV.write(snakemake.output[2], sum_stats)
else

    #mpaths = choose_optimal_codepaths(pnts, cb, H, params, cpaths_df)

    # Search for codepaths in each cell to reduce memory allocation
    #cell_pnts = DataFrame.(collect(groupby(pnts, :cellid)))

    cell_cpaths = DataFrame.(collect(groupby(cpaths_df, :cellid)))
    cell_tups = [(pnts[pnts.cellid .== cell_cpaths[i].cellid[1],:], cell_cpaths[i]) for i in 1:length(cell_cpaths)]

    function choose_cell_optimal_codepaths(cell_tup)
        cell_pnts = cell_tup[1]
        cell_cpaths = cell_tup[2]

        #cell_mpaths, dense_discarded_cpaths = choose_optimal_codepaths(cell_pnts, cb, H, params, cell_cpaths, HiGHS.Optimizer, ret_discarded=true)

        cell_mpaths, dense_discarded_cpaths = choose_optimal_codepaths(cell_pnts, cb, H, params, cell_cpaths, Gurobi.Optimizer, ret_discarded=true)
        if nrow(cell_mpaths) > 0
            cell_mpaths[!,"cellid"] .= cell_cpaths.cellid[1]
        end
        return cell_pnts, cell_mpaths, dense_discarded_cpaths
    end

    cell_out = map(choose_cell_optimal_codepaths, cell_tups)
    cell_pnts = [cell[1] for cell in cell_out]
    cell_mpaths = [cell[2] for cell in cell_out]
    cell_discared_cpaths = [cell[3] for cell in cell_out]

    decoded_points = vcat(cell_pnts...)
    mpaths = vcat(cell_mpaths...)
    discared_cpaths = vcat(cell_discared_cpaths...)

    if nrow(mpaths) == 0
        CSV.write(snakemake.output[1], DataFrame(Dict("cpath"=>[],"cost"=>[],"gene_number"=>[],"x"=>[],"y"=>[],"z"=>[],"cellid"=>[])))
        sum_stats = DataFrame()
        #sum_stats[!, "pos"] = [pos]
        #sum_stats[!, "ch"] = [ch]
        sum_stats[!,"lvf"] = [lat_var_factor]
        sum_stats[!,"lwvf"] = [lw_var_factor]
        sum_stats[!,"svf"] = [s_var_factor]
        sum_stats[!,"dr"] = [ndrops]
        sum_stats[!,"n_barcodes"] = [0]
        sum_stats[!,"negative_control_barcodes"] = [0]
        sum_stats[!,"gene_encoding_barcodes"] = [0]
        sum_stats[!,"est_fp_rate"] = [0]
        CSV.write(snakemake.output[2], sum_stats)
        exit()
    end

    #pos = parse(Int64, snakemake.wildcards["pos"])
    #ch = parse(Int64, snakemake.wildcards["ch"])

    #mpaths[!,"pos"] .= pos
    #mpaths[!,"ch"] .= ch
    mpaths[!,"lvf"] .= lat_var_factor
    mpaths[!, "wf"] .= lw_var_factor
    mpaths[!, "svf"] .= s_var_factor
    mpaths[!, "dr"] .= ndrops
    sort!(mpaths, ["cellid", "cc"])

    n_codewords = nrow(cb)
    n_ontargets = n_codewords - sum(cb.gene .== "negative_control")

    n_barcodes = nrow(mpaths)
    n_neg_cntrl_mpaths = sum(mpaths.gene .== "negative_control")

    sum_stats = DataFrame()
    #sum_stats[!, "pos"] = [pos]
    #sum_stats[!, "ch"] = [ch]
    sum_stats[!,"lvf"] = [lat_var_factor]
    sum_stats[!,"lwvf"] = [lw_var_factor]
    sum_stats[!,"svf"] = [s_var_factor]
    sum_stats[!,"dr"] = [ndrops]
    sum_stats[!,"n_barcodes"] = [n_barcodes]
    sum_stats[!,"negative_control_barcodes"] = [n_neg_cntrl_mpaths]
    sum_stats[!,"gene_encoding_barcodes"] = [n_barcodes - n_neg_cntrl_mpaths]
    sum_stats[!,"est_fp_rate"] = [(n_neg_cntrl_mpaths*n_ontargets/(n_codewords-n_ontargets))/(n_barcodes-n_neg_cntrl_mpaths)]

    CSV.write(snakemake.output[1], mpaths)
    CSV.write(snakemake.output[2], sum_stats)
    #CSV.write(snakemake.output[3], decoded_points)
    #CSV.write(snakemake.output[4], discared_cpaths)
end
