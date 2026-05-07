using DataFrames, Parquet, GLMNet, SparseArrays, CSV, Plots, FileIO, Images

fname_cpaths = snakemake.input[1] #"results/cand_cpaths_cell_1.csv"
fname_cdots = snakemake.input[2] #"results/cand_dots/cand_dots_cell_1.csv"
fnames_images = snakemake.input[3] #"results/cell_imgs/aligned_stack_bgsub_ch_488_pos_0_cell_1.tif"

σ = snakemake.params["sigma"]
psf_cutoff = snakemake.params["psf_cutoff"]
min_peak = snakemake.params["min_peak"]

psf(x,y,z) = exp(-(x^2 + y^2 + z^2)/(2.0* σ^2))
psf(x,y) = exp(-(x^2 + y^2)/(2.0* σ^2))
psf(xyz ::Tuple) = psf(xyz...)

psf_mat = psf.(Tuple.(collect(CartesianIndices((-3:3, -3:3)))))
psf_mat[psf_mat .< psf_cutoff] .= 0

psf_is = Int32[]
psf_js = Int32[]
psf_vs = Float64[]
for i in 1:7, j in 1:7
    if psf_mat[i,j] > psf_cutoff
        push!(psf_is, Int32(i-4))
        push!(psf_js, Int32(j-4))
        push!(psf_vs, psf_mat[i,j])
    end
end

sparse_psf_df = DataFrame(Dict("im_row" => psf_is, "im_col" => psf_js, "v" => psf_vs))



function build_model(_cpaths_df, cand_dots, image_dict, sparse_psf_df, npcs)
    images_dims = {hyb_key => size(image_dict[hyb_key]) for hyb_key in keys(image_dict)}
        
    #row_dim, col_dim, zdim = size(hyperstack)
    npix_per_img_dict = {hyb_key => prod(images_dims[hyb_key]) for hyb_key in keys(images_dims)}   #row_dim*col_dim

    sparse_dots = DataFrame[]
    for (j, cpath) in enumerate(eachrow(_cpaths_df))
        for dot in cpath["cpath"]
            # copy psf vals
            sparse_dot_rep = copy(sparse_psf_df)
            
            # move to dot location
            #dot_x, dot_y, block, pseudocolor = cand_dots[dot, ["x", "y", "hyb"]] #"block", "pseudocolor"]]
            dot_x, dot_y, hyb = cand_dots[dot, ["x", "y", "hyb"]] #"block", "pseudocolor"]]
            sparse_dot_rep[!, "im_row"] .+= dot_x
            sparse_dot_rep[!, "im_col"] .+= dot_y

            (row_dim, col_dim) = images_dims[hyb]
            img_npix = npix_per_img_dict[hyb]
            
            # remove out of bound pixels
            filter!(p -> (p.im_row > 0) && (p.im_row <= row_dim) && (p.im_col > 0) && (p.im_col <= col_dim), sparse_dot_rep)

            # convert row-col of sparse block-pseudocolor image represntation to i of sparse design matrix
            design_mat_row = sparse_dot_rep.im_row .+ row_dim*(sparse_dot_rep.im_col .- 1)
            design_mat_row .+=  img_npix * (hyb -1 ) #npcs * (block - 1) + (pseudocolor - 1))

            

            push!(sparse_dots, DataFrame(Dict("row" => design_mat_row, "col" => fill(j, length(design_mat_row)), "v" => sparse_dot_rep.v)))
        end
    end
    Adf = vcat(sparse_dots...)
    
    # reindex rows
    ys = vcat([rehape(image_dict[hyb_key], npix_per_img[hyb_key]) for hyb_key in keys(image_dict)]...)
    #ys = reshape(hyperstack, prod(size(hyperstack)))
    sorted_unique_rows = Int64.(sort(unique(Adf.row)))
    println(typeof(sorted_unique_rows))
    ydf = DataFrame(Dict("row" => sorted_unique_rows, "y" => ys[sorted_unique_rows]))
        
    reindex_row_dict = Dict([(row0 => rowp) for (rowp, row0) in enumerate(sorted_unique_rows)])
    reindex_row(row) = reindex_row_dict[row]
    Adf[!,"row"] .= reindex_row.(Adf[!,"row"])
    ydf[!,"row"] .= reindex_row.(ydf[!,"row"])

    filter!(r -> r.y > 0, ydf)
    ydf = sort(unique(ydf))
    return ydf, Adf
end


# read in data
cpaths_df = DataFrame(CSV.File(fname_cpaths))
cpaths_df[!,"cpath"] = eval.(Meta.parse.(cpaths_df[!,"cpath"]))

candidate_dot_coords = DataFrame(read_parquet(fname_cdots))


candidate_dot_coords[!, "x"] .= Float64.(candidate_dot_coords[!, "x"] )
candidate_dot_coords[!, "y"] .= Float64.(candidate_dot_coords[!, "y"] )
candidate_dot_coords[!, "z"] .= Float64.(candidate_dot_coords[!, "z"] )
candidate_dot_coords[!, "s"] .= Float64.(candidate_dot_coords[!, "s"] )

candidate_dot_coords[!, "s_z"] .= Float64.(candidate_dot_coords[!, "s_z"] )
candidate_dot_coords[!, "s_xy"] .= Float64.(candidate_dot_coords[!, "s_xy"] )
candidate_dot_coords[!, "w"] .= Float64.(candidate_dot_coords[!, "w"] )
candidate_dot_coords[!, "block"] .= Int64.(candidate_dot_coords[!, "block"] )
candidate_dot_coords[!, "pseudocolor"] .= Int64.(candidate_dot_coords[!, "pseudocolor"] )
candidate_dot_coords[candidate_dot_coords[!, "pseudocolor"] .== 0, "pseudocolor"] .= 20

images = {split(fname, "_")[4] => round.(UInt16, reinterpret.(Float64, channelview(load(fname)))) for (i, fname) in enumerate(snakemake.input["imgs"])} 


#hyperstack = reinterpret.(Float64, channelview(load(fname_hstack)))
#hyperstack = round.(UInt16, hyperstack)

println("buidling model")
ydf, Adf = build_model(cpaths_df, candidate_dot_coords, images, sparse_psf_df, 20)

Y = SparseVector(maximum(Adf.row), UInt32.(ydf.row), UInt16.(ydf.y))
A = sparse(UInt32.(Adf.row), UInt32.(Adf.col), Float64.(Adf.v))
println("done")

constraints = zeros(2, maximum(Adf.col))
constraints[2,:] .= 2^16

println("running Model")
lasso_path = glmnet(A, Y, constraints = constraints)
println("done")

summary_df = DataFrame(Dict("lambda" => lasso_path.lambda, "intercept" => lasso_path.a0, "R_square" => lasso_path.dev_ratio))


ΔRsq =  lasso_path.dev_ratio[2:end] - lasso_path.dev_ratio[1:end-1]
ΔRsq_max = argmax(ΔRsq)
thresh_ind = findfirst(r -> r < 0.001, ΔRsq[ΔRsq_max+1:end])
chosen_coeff_ind = isnothing(thresh_ind) ? size(lasso_path.betas)[2] : ΔRsq_max + thresh_ind -1

println("thresholding...")
thresh_coeffs = lasso_path.betas[:, chosen_coeff_ind] .> min_peak
betas = lasso_path.betas[thresh_coeffs, chosen_coeff_ind]

beta_path = lasso_path.betas[thresh_coeffs, :]

println("subsetting...")
cpaths_df_select = cpaths_df[collect(1:length(thresh_coeffs))[thresh_coeffs], :]

lamba_rank_most_stringent = map(s -> findfirst(b -> b > 0, s), eachslice(beta_path, dims=1))
cpaths_df_select[!, "lambda_most_stringent"] = lasso_path.lambda[lamba_rank_most_stringent]
#cpaths_df_select[!, "beta_least_stringent"] = betas


println("saving")
CSV.write(snakemake.output[1], cpaths_df_select)
CSV.write(snakemake.output[2], summary_df)
save(snakemake.output[3], round.(UInt16, beta_path))
println("done")
