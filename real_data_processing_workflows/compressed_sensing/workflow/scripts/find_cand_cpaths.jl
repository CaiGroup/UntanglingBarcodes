using DataFrames, CSV, SeqFISHSyndromeDecoding, DelimitedFiles, Parquet, Pkg
Pkg.activate(snakemake.input["env"])


lat_search_radius = snakemake.params["lat_search_radius"] #1.1
z_search_radius = snakemake.params["z_search_radius"] #1
zeros_probed = true

cb = DataFrame(CSV.File(snakemake.input[1])) #"codebook_ch_488.csv"))
cb[!, "block1"] = Int8.(cb[!, "block1"] .% 20)
cb[!, "block2"] = Int8.(cb[!, "block2"] .% 20)
cb[!, "block3"] = Int8.(cb[!, "block3"] .% 20)
cb[!, "block4"] = Int8.(cb[!, "block4"] .% 20)


H = readdlm(snakemake.input[2], '\t', Int8) #readdlm("H.txt", '\t', Int8)

cell_num = 1

#candidate_dot_coords = DataFrame(CSV.File(snakemake.input[3])) #"cand_dot_coords_cell$cell_num.csv"))
candidate_dot_coords = DataFrame(read_parquet(snakemake.input[3])) #"cand_dot_coords_cell$cell_num.csv"))


candidate_dot_coords[!, "x"] .= Float64.(candidate_dot_coords[!, "x"] )
candidate_dot_coords[!, "y"] .= Float64.(candidate_dot_coords[!, "y"] )
candidate_dot_coords[!, "z"] .= Float64.(candidate_dot_coords[!, "z"] )
candidate_dot_coords[!, "s"] .= Float64.(candidate_dot_coords[!, "s"] )

candidate_dot_coords[!, "s_z"] .= Float64.(candidate_dot_coords[!, "s_z"] )
candidate_dot_coords[!, "s_xy"] .= Float64.(candidate_dot_coords[!, "s_xy"] )
candidate_dot_coords[!, "w"] .= Float64.(candidate_dot_coords[!, "w"] )
candidate_dot_coords[!, "block"] .= Int64.(candidate_dot_coords[!, "block"] )
candidate_dot_coords[!, "pseudocolor"] .= Float64.(candidate_dot_coords[!, "pseudocolor"] )
candidate_dot_coords[!, "hyb"] .= Int64.(candidate_dot_coords[!, "hyb"] )


#candidate_dot_coords[!, "w"] = zeros(nrow(candidate_dot_coords))
#candidate_dot_coords[!, "s"] = zeros(nrow(candidate_dot_coords))


params = DecodeParams()
set_xy_search_radius(params, lat_search_radius)
set_z_search_radius(params, z_search_radius)
set_n_allowed_drops(params, 0)
set_lat_var_cost_coeff(params, 0)
set_z_var_cost_coeff(params, 0)
set_lw_var_cost_coeff(params, 0)
set_s_var_cost_coeff(params, 0)
set_zeros_probed(params, zeros_probed)

candidate_dot_coords[!, "dotID"] .= 1:nrow(candidate_dot_coords)
println("finding codepaths")
all_codepaths = get_codepaths(candidate_dot_coords, cb, H, params)
println("done: saving...")
CSV.write(snakemake.output[1], all_codepaths)
println("saved")
#CSV.write(snakemake.output[1], all_codepaths)

