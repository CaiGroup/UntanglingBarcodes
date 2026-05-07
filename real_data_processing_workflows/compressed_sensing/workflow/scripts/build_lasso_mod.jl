using Pkg
Pkg.activate(snakemake.input["env"])

using CSV, DataFrames, GLMNet, SeqFISHSyndromeDecoding, Parquet, SparseArrays, Images, StatsBase
println("snakemake.wildcards: ", snakemake.wildcards)

#fname_cpaths = snakemake.input["cand_cpaths"] #"results/cand_cpaths_cell_1.csv"
#fname_cdots = snakemake.input["cand_dots_unregistered"] #"results/cand_dots/cand_dots_cell_1.csv"
#fname_hstack = snakemake.input[3] #"results/cell_imgs/aligned_stack_bgsub_ch_488_pos_0_cell_1.tif"

#σ = snakemake.params["sigma"]
psf_cutoff = snakemake.params["psf_cutoff"]
min_peak = snakemake.params["min_peak"]
sr = parse(Float64, snakemake.wildcards["sr"])
#ch = parse(Int64, snakemake.wildcards["ch"])
ch = snakemake.wildcards["ch"]
σ = snakemake.params["sigma"][ch]

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
    images_dims = Dict((hyb_key => size(image_dict[hyb_key])) for hyb_key in keys(image_dict))
    sorted_hyb_keys = sort(collect(keys(images_dims)))
    #println("sorted_hyb_keys: ", sorted_hyb_keys)
        
    #row_dim, col_dim, zdim = size(hyperstack)
    npix_per_img_dict = Dict((hyb_key => prod(images_dims[hyb_key])) for hyb_key in keys(images_dims))   #row_dim*col_dim

    img_sizes = [npix_per_img_dict[hyb_key] for hyb_key in sorted_hyb_keys]
    #println("get cum sums")
    #println(typeof(img_sizes))
    cum_img_sizes = cumsum(img_sizes)
    insert!(cum_img_sizes, 1, 0)
    #println("build design matrix")

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
            #img_npix = npix_per_img_dict[hyb]
            
            # remove out of bound pixels
            filter!(p -> (p.im_row > 0) && (p.im_row <= row_dim) && (p.im_col > 0) && (p.im_col <= col_dim), sparse_dot_rep)

            # convert row-col of sparse block-pseudocolor image represntation to i of sparse design matrix
            design_mat_row = sparse_dot_rep.im_row .+ row_dim*(sparse_dot_rep.im_col .- 1)
            design_mat_row .+= cum_img_sizes[hyb+1] #npcs * (block - 1) + (pseudocolor - 1))

            

            push!(sparse_dots, DataFrame(Dict("row" => design_mat_row, "col" => fill(j, length(design_mat_row)), "v" => sparse_dot_rep.v)))
        end
    end
    Adf = vcat(sparse_dots...)
    
    # reindex rows
    ys = vcat([reshape(image_dict[hyb_key], npix_per_img_dict[hyb_key]) for hyb_key in sort(collect(keys(image_dict)))]...)
    sorted_unique_rows = Int64.(sort(unique(Adf.row)))
    ydf = DataFrame(Dict("row" => sorted_unique_rows, "y" => Int64.(ys[sorted_unique_rows])))
    reindex_row_dict = Dict([(row0 => rowp) for (rowp, row0) in enumerate(sorted_unique_rows)])
    reindex_row(row) = reindex_row_dict[row]
    Adf[!,"row"] .= reindex_row.(Adf[!,"row"])
    ydf[!,"row"] .= reindex_row.(ydf[!,"row"])
    filter!(r -> r.y > 0, ydf)
    ydf = sort(unique(ydf))
    return ydf, Adf
end

#params = DecodeParams()

#set_xy_search_radius(params, 4.0)
#set_z_search_radius(params, 1.0)
#set_n_allowed_drops(params, 0)
#set_lat_var_cost_coeff(params, 1.0)
#set_z_var_cost_coeff(params, 0.0)
#set_lw_var_cost_coeff(params, 0.0)
#set_s_var_cost_coeff(params, 0.0)
#set_zeros_probed(params, true)


# read in registered data (dots aligned with each other)
candidate_dot_coords = DataFrame(read_parquet(snakemake.input["cand_dots_aligned"]))
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

#cost(cpath) = SeqFISHSyndromeDecoding.obj_function(cpath, candidate_dot_coords, 4, params, nothing)

cpaths_df = DataFrame(CSV.File(snakemake.input["cand_cpaths"]))
cpaths_df[!,"cpath"] = eval.(Meta.parse.(cpaths_df[!,"cpath"]))
cpaths_df = SeqFISHSyndromeDecoding.threshold_cpaths(cpaths_df, candidate_dot_coords, sr, 1.0, nothing)


# read in unregistered data (dots aligned with images)
candidate_dot_coords = DataFrame(read_parquet(snakemake.input["cand_dots_unregistered"]))
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

#cost(cpath) = SeqFISHSyndromeDecoding.obj_function(cpath, candidate_dot_coords, 4, params, nothing)
#cpaths_df[!, "cost"] = lat_var_penalty .* cost.(cpaths_df[!, "cpath"]) .+ 1

#images = {split(fname, "_")[4] => round.(UInt16, reinterpret.(Float64, channelview(load(fname)))) for (i, fname) in enumerate(snakemake.input["imgs"])} 
images = Dict((parse(Int64, split(fname, "_")[4]) => reinterpret.(UInt16, channelview(load(fname)))) for (i, fname) in enumerate(snakemake.input["imgs"]))


#hyperstack = reinterpret.(Float64, channelview(load(fname_hstack)))
#hyperstack = round.(UInt16, hyperstack)

cpaths_no_neg_ctrl = filter(cp -> cp.gene != "negative_control", cpaths_df)

ydf_wnc, Adf_wnc = build_model(cpaths_df, candidate_dot_coords, images, sparse_psf_df, 20)


cpaths_no_neg_ctrl = filter(cp -> cp.gene != "negative_control", cpaths_df)

ydf_nnc, Adf_nnc = build_model(cpaths_no_neg_ctrl, candidate_dot_coords, images, sparse_psf_df, 20)


write_parquet(snakemake.output[1], ydf_wnc)
write_parquet(snakemake.output[2], Adf_wnc)
CSV.write(snakemake.output[3], cpaths_df)
write_parquet(snakemake.output[4], ydf_nnc)
write_parquet(snakemake.output[5], Adf_nnc)
CSV.write(snakemake.output[6], cpaths_no_neg_ctrl)
