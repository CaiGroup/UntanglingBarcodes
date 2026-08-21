using Pkg
Pkg.activate(snakemake.input[4])

using SeqFISHSyndromeDecoding
using CSV
using DataFrames
using Statistics
using Profile
using DelimitedFiles
using Statistics

code = snakemake.wildcards["code"]

hyb_pnts_filename = snakemake.input[1]#"alignment/aligned_dots_mw_600_ch_0.csv"

println("read pnts")
pnts = DataFrame(CSV.File(hyb_pnts_filename))#
pnts.block = UInt8.(pnts.block)

#if typeof(pnts.pseudocolor[1]) != String7 || typeof(pnts.pseudocolor[1]) != String
pnts.pseudocolor = UInt8.(pnts.pseudocolor)
pnts.x = Float64.(pnts.x)
pnts.y = Float64.(pnts.y)
pnts.z = zeros(Float64, nrow(pnts))

#H = readdlm(snakemake.input[3], ',', Int64)
H = readdlm(snakemake.input[3], ',')

if any(typeof.(H) .<: AbstractString) # == String7 || typeof(H[1,1]) == String || typeof(H[1,1]) == SubString{String}
    H = string.(H)
elseif any(H .< 0)
    H = Int8.(H)
else
    H = UInt8.(H)
end

#H = [1 1 1 -1;]


#lat_thresh = 2.0
z_thresh = 0.0
#code = "E2019"#_neg_ctrl"
#code = "E2019_neg_ctrl"
cb_name = snakemake.input[2] #"E2019_cb_all_control.txt"
#cb_name = "E2019_cb_all_control.txt"
#cb = readdlm("codebooks/" * cb_name, UInt8)
#cb = readdlm(cb_name, UInt8)
cb = DataFrame(CSV.File(cb_name))

#ndrops = 0#parse(Int64, snakemake.wildcards["dr"])#1
ndrops = snakemake.params["drc"]
#ndrops = 0

println("add dot id")
pnts.dot_ID = Array(1:nrow(pnts))
ndots = nrow(pnts)

# Cost Parameters
free_dot_cost = 1.0
lat_var_factor = snakemake.params["lf"]
#lat_var_factor = 40.0#parse(Float64, snakemake.wildcards["lf"])
z_var_factor = 0.0
lw_var_factor = 0.0 #snakemake.params["wf"]
#lw_var_factor = 8.0#parse(Float64, snakemake.wildcards["wf"])
s_var_factor = 0.0 #snakemake.params["sf"]
#s_var_factor = 4.0#parse(Float64, snakemake.wildcards["sf"])

erasure_penalty = 2*free_dot_cost

#lat_thresh = sqrt(20.0/lat_var_factor)*4
#lat_thresh = sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*4 #used for 20210528_results
lat_thresh = snakemake.params["lat_thresh"] #sqrt(free_dot_cost*size(H)[2]/lat_var_factor)*3
rstdv = parse(Float64, snakemake.wildcards["rstdv"])


println("cost parameters: $lat_var_factor, $lw_var_factor")

params = DecodeParams()
set_free_dot_cost(params, free_dot_cost)
set_xy_search_radius(params, rstdv * lat_thresh)
set_z_search_radius(params, z_thresh)
set_n_allowed_drops(params, ndrops)
set_lat_var_cost_coeff(params, lat_var_factor)
set_z_var_cost_coeff(params, z_var_factor)
set_lw_var_cost_coeff(params, lw_var_factor)
set_s_var_cost_coeff(params, s_var_factor)
set_erasure_penalty(params, erasure_penalty) #snakemake.params["drop_penalty"])
print("code $code")
if code[1:7] == "seqFISH"
    println("Zeros Probed")
    set_zeros_probed(params, true)
else
    println("zeros not probed")
    set_zeros_probed(params, false)
end
set_skip_thresh(params, 10000000)
set_skip_density_thresh(params, 100000000000)

# Search for codepaths in each cell to reduce memory allocation
function get_cell_codepaths(pnts_)
    cpaths = get_codepaths(pnts_, cb, H, params)
    ncpaths = typeof(cpaths) == DataFrame ? nrow(cpaths) : length(cpaths)
    if ncpaths == 0
	return
    end
    insertcols!(cpaths, :cellid=>fill(pnts_.cellid[1], ncpaths))
    return cpaths
end

cell_pnts = DataFrame.(collect(groupby(pnts, :cellid)))
println("got cell points")
filter!(df -> nrow(df) >= (4-ndrops), cell_pnts)
cell_cpaths = map(get_cell_codepaths, cell_pnts)

filter!(df -> typeof(df) == DataFrame, cell_cpaths)
cpaths = vcat(cell_cpaths...)
#pnts2 = vcat(cell_pnts...)
println("saving...")
if typeof(cpaths) != DataFrame
    CSV.write(snakemake.output[1], DataFrame(Dict("cpath"=>[],"cost"=>[],"gene_number"=>[],"x"=>[],"y"=>[],"z"=>[],"cellid"=>[])))
elseif nrow(cpaths) > 2^22
    CSV.write(snakemake.output[1],cpaths,bufsize=nrow(cpaths)+1000)
else
    CSV.write(snakemake.output[1], cpaths)
end
