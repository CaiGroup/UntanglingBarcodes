using Pkg
Pkg.activate(snakemake.input["env"])
using CSV, DataFrames, GLMNet, SeqFISHSyndromeDecoding, Parquet, SparseArrays, Images, StatsBase

σ = snakemake.params["sigma"]
psf_cutoff = snakemake.params["psf_cutoff"]
#ch = parse(Int64, snakemake.wildcards["ch"])
ch = snakemake.wildcards["ch"]
min_peak = snakemake.params["min_peak"][ch]
sr = parse(Float64, snakemake.wildcards["sr"])
lat_var_penalty = parse(Float64, snakemake.wildcards["lf"])

function run_lasso(_cpaths_df, ydf, Adf, min_peak; intercept=true, lambdas = nothing)

    Y = SparseVector(maximum(Adf.row), UInt32.(ydf.row), UInt16.(ydf.y))
    A = sparse(UInt32.(Adf.row), UInt32.(Adf.col), Float64.(Adf.v))

    constraints = zeros(2, maximum(Adf.col))
    constraints[2,:] .= 2^16

    if isnothing(lambdas)
        lasso_path = glmnet(A, Y, constraints = constraints, penalty_factor = _cpaths_df[!, "cost"]; intercept=intercept)
    else
        lasso_path = glmnet(A, Y, constraints = constraints, penalty_factor = _cpaths_df[!, "cost"]; intercept=intercept, lambda=lambdas)
    end

    _cpaths_df[!, "largest_lambda_over_thresh"] = map(s -> findfirst(b -> b > min_peak, s), eachslice(lasso_path.betas, dims=1))
    beta_paths_selected = lasso_path.betas[ .~ isnothing.(_cpaths_df[!, "largest_lambda_over_thresh"]), :]
    selected_cps = _cpaths_df[ .~ isnothing.(_cpaths_df[!, "largest_lambda_over_thresh"]), :]
    selected_cps[!, "largest_lambda_over_thresh"] = lasso_path.lambda[selected_cps[!, "largest_lambda_over_thresh"]]


    lambda_marge_nonzeros = countmap(selected_cps[!, "largest_lambda_over_thresh"])
    #df = cumsum(map(p -> p.second, sort(collect(lambda_marge_nonzeros), rev=true)))
    lasso_path_summary = DataFrame(Dict("lambda" => lasso_path.lambda, "intercept" => lasso_path.a0, "R_square" => lasso_path.dev_ratio)) #, "df" =>df ))



    return selected_cps, lasso_path_summary, beta_paths_selected
end

params = DecodeParams()

set_xy_search_radius(params, sr)
set_z_search_radius(params, 1.0)
set_n_allowed_drops(params, 0)
set_lat_var_cost_coeff(params, 1.0)
set_z_var_cost_coeff(params, 0.0)
set_lw_var_cost_coeff(params, 0.0)
set_s_var_cost_coeff(params, 0.0)
set_erasure_penalty(params, 1000.0)
set_zeros_probed(params, true)


# read in registered data (dots aligned with each other)
candidate_dot_coords = DataFrame(read_parquet(snakemake.input["cdots_aligned"]))
candidate_dot_coords[!, "x"] .= Float64.(candidate_dot_coords[!, "x"] )
candidate_dot_coords[!, "y"] .= Float64.(candidate_dot_coords[!, "y"] )
candidate_dot_coords[!, "z"] .= Float64.(candidate_dot_coords[!, "z"] )
candidate_dot_coords[!, "s"] .= Float64.(candidate_dot_coords[!, "s"] )

candidate_dot_coords[!, "s_z"] .= Float64.(candidate_dot_coords[!, "s_z"] )
candidate_dot_coords[!, "s_xy"] .= Float64.(candidate_dot_coords[!, "s_xy"] )
candidate_dot_coords[!, "w"] .= Float64.(candidate_dot_coords[!, "w"] )
candidate_dot_coords[!, "block"] .= Int64.(candidate_dot_coords[!, "block"] )
candidate_dot_coords[!, "pseudocolor"] .= Int64.(candidate_dot_coords[!, "pseudocolor"] )
candidate_dot_coords[!, "hyb"] .= Int64.(candidate_dot_coords[!, "hyb"] )
candidate_dot_coords[candidate_dot_coords[!, "pseudocolor"] .== 0, "pseudocolor"] .= 20

cost(cpath) = SeqFISHSyndromeDecoding.obj_function(cpath, candidate_dot_coords, 4, params, nothing)

cpaths_df_nnc = DataFrame(CSV.File(snakemake.input["cpaths_nnc"]))
cpaths_df_nnc[!,"cpath"] = eval.(Meta.parse.(cpaths_df_nnc[!,"cpath"]))
#println("pre threshold nrow(cpaths_df) = ", nrow(cpaths_df))
#cpaths_df = SeqFISHSyndromeDecoding.threshold_cpaths(cpaths_df, candidate_dot_coords, sr, 1.0, nothing)
#println("post threshold nrow(cpaths_df) = ", nrow(cpaths_df))

cpaths_df_nnc[!, "cost"] = lat_var_penalty .* cost.(cpaths_df_nnc[!, "cpath"]) .+ 1


Y_nnc = read_parquet(snakemake.input["y_nnc"])
A_nnc = read_parquet(snakemake.input["A_nnc"])

# solve for case of no negative controls
cp_select_nnc, pth_sm_nnc, betas_nnc = run_lasso(cpaths_df_nnc, Y_nnc, A_nnc, min_peak)


cpaths_df_wnc = DataFrame(CSV.File(snakemake.input["cpaths_wnc"]))
cpaths_df_wnc[!,"cpath"] = eval.(Meta.parse.(cpaths_df_wnc[!,"cpath"]))
Y_wnc = read_parquet(snakemake.input["y_wnc"])
A_wnc = read_parquet(snakemake.input["A_wnc"])

cpaths_df_wnc[!, "cost"] = lat_var_penalty .* cost.(cpaths_df_wnc[!, "cpath"]) .+ 1

# evaulate with negative controls at the same lambdas as without negative controls
cp_select_wnc, pth_sm_wnc, betas_wnc = run_lasso(cpaths_df_wnc, Y_wnc, A_wnc, min_peak, lambdas=pth_sm_nnc.lambda)

rename!(pth_sm_nnc, :R_square => :R_square_nnc)
rename!(pth_sm_nnc, :intercept => :intercept_nnc)
#rename!(pth_sm_nnc, :df => :df_nnc)

pth_sm_nnc[!, "R_square_wnc"] = pth_sm_wnc[!, "R_square"]
pth_sm_nnc[!, "intercept_wnc"] = pth_sm_wnc[!, "intercept"]
#pth_sm_nnc[!, "df_wnc"] = pth_sm_wnc[!, "df"]

CSV.write(snakemake.output[1], cp_select_nnc)
CSV.write(snakemake.output[2], cp_select_wnc)
CSV.write(snakemake.output[3], pth_sm_nnc)
save(snakemake.output[4], round.(UInt16, betas_nnc))
save(snakemake.output[5], round.(UInt16, betas_wnc))
