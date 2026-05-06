using Pkg
#Pkg.activate(snakemake.input[4])
using CSV, DataFrames, DelimitedFiles, SeqFISHSyndromeDecoding

cb = DataFrame(CSV.File(snakemake.input[1]))
H = readdlm(snakemake.input[2], ',', UInt8)

candidate_dot_coords = DataFrame(CSV.File(snakemake.input[3]))
SeqFISHSyndromeDecoding.sort_readouts!(candidate_dot_coords)
lat_search_radius = snakemake.params["sr"]
z_search_radius = 1
zeros_probed = false

function fcpths(candidate_dot_coords, cb, H, lat_search_radius, z_search_radius, zeros_probed)
    params = DecodeParams()
    set_xy_search_radius(params, lat_search_radius)
    set_z_search_radius(params, z_search_radius)
    set_n_allowed_drops(params, 0)
    set_lat_var_cost_coeff(params, 0)
    set_z_var_cost_coeff(params, 0)
    set_lw_var_cost_coeff(params, 0)
    set_s_var_cost_coeff(params, 0)
    set_zeros_probed(params, zeros_probed)

    codepaths = get_codepaths(candidate_dot_coords, cb, H, params)
    
    return codepaths
end

cpaths = fcpths(candidate_dot_coords, cb, H, lat_search_radius, z_search_radius, zeros_probed)

CSV.write(snakemake.output[1], cpaths)

