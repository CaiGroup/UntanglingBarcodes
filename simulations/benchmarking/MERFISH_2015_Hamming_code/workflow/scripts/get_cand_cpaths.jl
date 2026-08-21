using Pkg
Pkg.activate(snakemake.input[4])
using CSV, DataFrames, DelimitedFiles, SeqFISHSyndromeDecoding

cb = DataFrame(CSV.File(snakemake.input[1]))
H = readdlm(snakemake.input[2], ',')
if any(typeof.(H) .<: AbstractString) # == String7 || typeof(H[1,1]) == String || typeof(H[1,1]) == SubString{String}
    H = string.(H)
elseif any(H .< 0)
    H = Int8.(H)
else
    H = UInt8.(H)
end

candidate_dot_coords = DataFrame(CSV.File(snakemake.input[3]))
SeqFISHSyndromeDecoding.sort_readouts!(candidate_dot_coords)
lat_search_radius = snakemake.params["sr"]
rstdv = parse(Float64, snakemake.wildcards["rstdv"])
z_search_radius = 1

code = snakemake.wildcards["code"]
if code[1:7] == "seqFISH"
    println("Zeros Probed")
    zeros_probed = true
else
    println("zeros not probed")
    zeros_probed = false
end
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

cpaths = fcpths(candidate_dot_coords, cb, H, rstdv * lat_search_radius, z_search_radius, zeros_probed)

CSV.write(snakemake.output[1], cpaths)
