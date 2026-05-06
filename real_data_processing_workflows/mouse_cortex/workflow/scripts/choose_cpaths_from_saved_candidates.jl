using Pkg
Pkg.activate(snakemake.input[5])

#using SeqFISHPointDecoding
using SeqFISHSyndromeDecoding
using CSV
using DataFrames
using Statistics
using Profile
using DelimitedFiles
using Gurobi
using FileIO
using Images
using GLPK

cpaths_fname = snakemake.input[1]
hyb_pnts_filename =  snakemake.input[2]#"alignment/aligned_dots_mw_600_ch_0.csv"
cb_name = snakemake.input[3] #"E2019_cb_all_control.txt"

pnts = DataFrame(CSV.File(hyb_pnts_filename))#
filter!(pnt -> ~ismissing(pnt.pseudocolor), pnts)
pnts.block = UInt8.(pnts.block)
pnts.x = Float64.(pnts.x)
pnts.y = Float64.(pnts.y)
#pnts.z = Float64.(pnts.z) #zeros(Float64, nrow(pnts))
pnts.z = zeros(Float64, nrow(pnts))
#pnts.s = pnts.sxy

select!(pnts, Not(:hyb))

select!(pnts, :block, :pseudocolor, :x, :y, :z, :s, :w, :cellid)
SeqFISHSyndromeDecoding.sort_readouts!(pnts)

function write_null()
    mpaths = DataFrame()
    
    mpaths[!, "gene"] = []
    mpaths[!, "gene_number"] = []
    mpaths[!, "cpath"] = []
    mpaths[!, "cost"] = []
    mpaths[!, "x"] = []
    mpaths[!, "y"] = []
    mpaths[!, "z"] = []
    mpaths[!, "cellid"] = []
    mpaths[!, "cc"] = []
    mpaths[!, "cc_size"] = []
    mpaths[!,"pos"] = []
    mpaths[!,"cell"] = []
    mpaths[!,"ch"] = []
    mpaths[!,"lvf"] = []
    mpaths[!,"zvf"] = []
    mpaths[!, "wf"] = []
    mpaths[!, "svf"] = []
    mpaths[!, "dr"] = []

    sum_stats = DataFrame()
    sum_stats[!, "pos"] = []
    sum_stats[!, "cell"] = []
    sum_stats[!, "ch"] = []
    #sum_stats[!, "mw"] = [mw]
    sum_stats[!,"lvf"] = []
    sum_stats[!,"zvf"] =[]
    sum_stats[!,"lwvf"] = []
    sum_stats[!,"svf"] = []
    sum_stats[!,"dr"] = []
    sum_stats[!,"n_barcodes"] = []
    sum_stats[!,"negative_control_barcodes"] = []
    sum_stats[!,"gene_encoding_barcodes"] = []
    sum_stats[!,"est_fp_rate"] = []

    CSV.write(snakemake.output[1], mpaths)
    CSV.write(snakemake.output[2], sum_stats)
    CSV.write(snakemake.output[4], DataFrame("na"=>0))
end

if nrow(pnts) == 0
    write_null()
    CSV.write(snakemake.output[3], DataFrame("na"=>0))
    exit()
else
    H = readdlm(snakemake.input[4], ',', UInt8)
    println("H")
    println(H)


    #z_thresh = 0.0
    cb = DataFrame(CSV.File(cb_name))

    println(first(cb,5))

    ndrops = parse(Int64, snakemake.wildcards["dr"])#1
    #mw = parse(Float64, snakemake.wildcards["mw"])


    pnts.dot_ID = Array(1:nrow(pnts))
    ndots = nrow(pnts)

    # Cost Parameters
    free_dot_cost = 1.0
    lat_var_factor = parse(Float64, snakemake.wildcards["lf"])

    z_var_factor = parse(Float64, snakemake.wildcards["zf"])
    lw_var_factor = parse(Float64, snakemake.wildcards["wf"])
    s_var_factor = parse(Float64, snakemake.wildcards["sf"])

    erasure_penalty = 4.0
    lat_thresh = snakemake.params["rxy"] #sqrt(free_dot_cost*4/lat_var_factor)*2
    z_thresh = snakemake.params["rz"] #sqrt(free_dot_cost*4/z_var_factor)*2


    println("cost parameters: $lat_var_factor, $lw_var_factor")

    params = DecodeParams()

    set_xy_search_radius(params, lat_thresh)
    set_z_search_radius(params, z_thresh)
    set_n_allowed_drops(params, ndrops)
    set_lat_var_cost_coeff(params, lat_var_factor)
    set_z_var_cost_coeff(params, z_var_factor)
    set_lw_var_cost_coeff(params, lw_var_factor)
    set_s_var_cost_coeff(params, s_var_factor)
    set_zeros_probed(params, false)
    if snakemake.params["skip_thresh"] == "auto"
        set_skip_thresh(params, nrow(cb))
    else
        set_skip_thresh(params, snakemake.params["skip_thresh"])
    end
    set_skip_density_thresh(params, snakemake.params["skip_density_thresh"])

    set_erasure_penalty(params, 2.0) #params, snakemake.params["drop_penalty"])

    cpaths_df = DataFrame(CSV.File(cpaths_fname))
    cpaths_df[!,"cpath"] = eval.(Meta.parse.(cpaths_df[!,"cpath"]))


    #mpaths = choose_optimal_codepaths(pnts, cb, H, params, cpaths_df)

    # Search for codepaths in each cell to reduce memory allocation
    #cell_pnts = DataFrame.(collect(groupby(pnts, :cellid)))


    cell_cpaths = DataFrame.(collect(groupby(cpaths_df, :cellid)))
    cell_tups = [(pnts[pnts.cellid .== cell_cpaths[i].cellid[1],:], cell_cpaths[i]) for i in 1:length(cell_cpaths)]

    function choose_cell_optimal_codepaths(cell_tup)
        cell_pnts = cell_tup[1]
        cell_cpaths = cell_tup[2]
        cell_mpaths, dense_discarded_cpaths = choose_optimal_codepaths(cell_pnts, cb, H, params, cell_cpaths, Gurobi.Optimizer, ret_discarded=true)
       if typeof(cell_mpaths) == DataFrame
            if nrow(cell_mpaths) > 0
                cell_mpaths[!,"cellid"] .= cell_cpaths.cellid[1]
            else
                cell_mpaths[!,"cellid"] .= [] #cell_cpaths.cellid[1]
            end
        end
        return cell_pnts, cell_mpaths, dense_discarded_cpaths
    end

    cell_out = map(choose_cell_optimal_codepaths, cell_tups)
    cell_pnts = [cell[1] for cell in cell_out]
    cell_mpaths = [cell[2] for cell in cell_out]
    cell_discarded_cpaths = [cell[3] for cell in cell_out]

    filter!(df -> typeof(df) == DataFrame, cell_pnts)
    filter!(df -> typeof(df) == DataFrame, cell_mpaths)
    filter!(df -> typeof(df) == DataFrame, cell_discarded_cpaths)

    decoded_points = vcat(cell_pnts...)
    mpaths = vcat(cell_mpaths...)
    discarded_cpaths = vcat(cell_discarded_cpaths...)


    #mpaths, discarded_cpaths = choose_optimal_codepaths(pnts, cb, H, params, cpaths_df, ret_discarded=true)
    #decoded_points = pnts

    pos = parse(Int64, snakemake.wildcards["pos"])
    cell = parse(Int64, snakemake.wildcards["cell"])
    #ch = parse(Int64, snakemake.wildcards["ch"])


    if typeof(mpaths) == DataFrame && nrow(mpaths) > 0
        println("nrow(mpaths): ", nrow(mpaths))
        mpaths[!,"pos"] .= pos
        mpaths[!,"cell"] .= cell 
        #mpaths[!,"ch"] .= ch
        mpaths[!,"lvf"] .= lat_var_factor
        mpaths[!,"zvf"] .= z_var_factor
        mpaths[!, "wf"] .= lw_var_factor
        mpaths[!, "svf"] .= s_var_factor
        mpaths[!, "dr"] .= ndrops
        #sort!(mpaths, ["cellid", "cc"])
        sort!(mpaths, "cc")

        n_codewords = nrow(cb)
        n_ontargets = n_codewords - sum(cb.gene .== "negative_control")

        n_barcodes = nrow(mpaths)
        n_neg_cntrl_mpaths = sum(mpaths.gene .== "negative_control")

        sum_stats = DataFrame()
        sum_stats[!, "pos"] = [pos]
        sum_stats[!, "cell"] = [cell]
        #sum_stats[!, "ch"] = [ch]
        #sum_stats[!, "mw"] = [mw]
        sum_stats[!,"lvf"] = [lat_var_factor]
        sum_stats[!,"zvf"] =[z_var_factor]
        sum_stats[!,"lwvf"] = [lw_var_factor]
        sum_stats[!,"svf"] = [s_var_factor]
        sum_stats[!,"dr"] = [ndrops]
        sum_stats[!,"n_barcodes"] = [n_barcodes]
        sum_stats[!,"negative_control_barcodes"] = [n_neg_cntrl_mpaths]
        sum_stats[!,"gene_encoding_barcodes"] = [n_barcodes - n_neg_cntrl_mpaths]
        sum_stats[!,"est_fp_rate"] = [(n_neg_cntrl_mpaths*n_ontargets/(n_codewords-n_ontargets))/(n_barcodes-n_neg_cntrl_mpaths)]

        CSV.write(snakemake.output[1], mpaths)
        CSV.write(snakemake.output[2], sum_stats)
        CSV.write(snakemake.output[3], decoded_points)
        CSV.write(snakemake.output[4], discarded_cpaths)
    else
        println("mpaths not dataframe:")
        println(mpaths)
        #CSV.write(snakemake.output[3], decoded_points)
        CSV.write(snakemake.output[3], DataFrame("na"=> 0))
        write_null()
    end
end

println("wrote decoded barcodes")
