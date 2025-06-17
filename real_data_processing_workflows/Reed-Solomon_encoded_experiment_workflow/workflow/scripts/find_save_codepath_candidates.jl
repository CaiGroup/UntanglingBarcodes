using SeqFISHSyndromeDecoding
using CSV
using DataFrames
using Statistics
using Profile
using DelimitedFiles
using Statistics
using Images
using FileIO

hyb_pnts_filename = snakemake.input[1]#"alignment/aligned_dots_mw_600_ch_0.csv"

pnts = DataFrame(CSV.File(hyb_pnts_filename))#
filter!(pnt -> ~ismissing(pnt.pseudocolor), pnts)
pnts.block = UInt8.(pnts.block)
pnts.x = Float64.(pnts.x)
pnts.y = Float64.(pnts.y)
#pnts.z = Float64.(pnts.z) #zeros(Float64, nrow(pnts))
pnts.z = zeros(Float64, nrow(pnts))
#pnts.s = pnts.sxy

select!(pnts, Not([:ch,:hyb]))
SeqFISHSyndromeDecoding.sort_readouts!(pnts)

if nrow(pnts) == 0
    CSV.write(snakemake.output[1], DataFrame(Dict("cpath"=>[],"cost"=>[],"gene_number"=>[],"x"=>[],"y"=>[],"z"=>[],"cellid"=>[])))
else

    filter!(pnt -> pnt.cellid != 0, pnts)
    #filter!(pnt -> pnt.cellid == 9, pnts)
    #filter!(pnt -> pnt.y > 850 && pnt.x > 850 && pnt.y < 900 && pnt.x < 900, pnts)
    println("npnts: ", nrow(pnts))

    H = readdlm(snakemake.input[3], ',', UInt8)

    #lat_thresh = 2.0
    #z_thresh = 0.0
    #code = "E2019"#_neg_ctrl"
    #code = "E2019_neg_ctrl"
    cb_name = snakemake.input[2] #"E2019_cb_all_control.txt"
    #cb_name = "E2019_cb_all_control.txt"!
    #cb = readdlm("codebooks/" * cb_name, UInt8)
    #cb = readdlm(cb_name, UInt8)
    cb = DataFrame(CSV.File(cb_name))

    #ndrops = 0#parse(Int64, snakemake.wildcards["dr"])#1
    ndrops = snakemake.params["dr"]
    println("ndrops: ", ndrops)
    #ndrops = 0


    pnts.dot_ID = Array(1:nrow(pnts))
    ndots = nrow(pnts)

    # Cost Parameters
    free_dot_cost = 1.0
    lat_var_factor = snakemake.params["lf"]
    ##lat_var_factor = 40.0#parse(Float64, snakemake.wildcards["lf"])
    z_var_factor = snakemake.params["zf"]
    #z_var_factor = parse(Float64,)

    lw_var_factor = snakemake.params["wf"]
    #lw_var_factor = 8.0#parse(Float64, snakemake.wildcards["wf"])
    s_var_factor = snakemake.params["sf"]
    #s_var_factor = 4.0#parse(Float64, snakemake.wildcards["sf"])

    erasure_penalty = 0.0

    #lat_thresh = sqrt(20.0/lat_var_factor)*4
    #lat_thresh = sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*4 #used for 20210528_results
    lat_thresh = snakemake.params["rxy"] #sqrt(free_dot_cost*4/lat_var_factor)*2 #sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*3
    z_thresh = snakemake.params["rz"] #sqrt(free_dot_cost*4/z_var_factor)*2 #sqrt(free_dot_cost*size(H)[2]/z_var_factor)*3


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
    #set_erasure_penalty(params, snakemake.params["drop_penalty"])

    # Search for codepaths in each cell to reduce memory allocation

    function get_cell_codepaths(pnts_)
        println("get codepaths of cell ", pnts_.cellid[1])
        cpaths = get_codepaths(pnts_, cb, H, params)
        if typeof(cpaths) == DataFrame && nrow(cpaths) > 1
            cpaths[!,"cellid"] .= pnts_.cellid[1]
            #cpaths[!,"x"] = mean.([pnts_.x[cpath] for cpath in cpaths.cpath])
            #cpaths[!,"y"] = mean.([pnts_.y[cpath] for cpath in cpaths.cpath])
            return cpaths
        end
    end

    println(first(pnts, 5))
    cell_pnts = DataFrame.(collect(groupby(pnts, :cellid)))
    println("length(cell_pnts): ", length(cell_pnts))
    cell_cpaths = map(get_cell_codepaths, cell_pnts)
    filter!(df -> typeof(df) == DataFrame, cell_cpaths)
    println("length filtered cpaths: ", length(cell_cpaths))

    #cpaths = get_codepaths(pnts, cb, H, params)
    cpaths = vcat(cell_cpaths...)
    #pnts2 = vcat(cell_pnts...)
    println("saving...")
    if typeof(cpaths) == DataFrame
        CSV.write(snakemake.output[1], cpaths)
    else
        CSV.write(snakemake.output[1], DataFrame(Dict("cpath"=>[],"cost"=>[],"gene_number"=>[],"x"=>[],"y"=>[],"z"=>[],"cellid"=>[])))
    end
    println("saved")
end
