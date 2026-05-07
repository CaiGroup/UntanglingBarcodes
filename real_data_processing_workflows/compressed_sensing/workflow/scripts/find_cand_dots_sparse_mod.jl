using Pkg
Pkg.activate(snakemake.input["env"])
using DataFrames, CSV, FileIO, Images, SparseArrays, GLMNet, Parquet

ch = parse(Int64, snakemake.wildcards["ch"])
ch = snakemake.wildcards["ch"]

println("ch: ", ch)
println(snakemake.params["sigma"])
σ = snakemake.params["sigma"][ch]
#σ = snakemake.params["sigma"]

hyb = parse(Int, snakemake.wildcards["hyb"])
pc = (hyb + 1) % 20 # parse(Int, snakemake.wildcards["pc"])
r = ceil((hyb + 1) / 20) #parse(Int, snakemake.wildcards["r"])
psf_cutoff = snakemake.params["psf_cutoff"]
min_peak = parse(Float64, snakemake.wildcards["mpcd"]) #snakemake.params["min_peak"][ch]
println("min_peak: ", min_peak)
psf(x,y,z) = exp(-(x^2 + y^2 + z^2)/(2.0* σ^2))

psf(x, y) = psf(x,y,0)
psf(t :: Tuple) = psf(t...)


get_conv_threshold(psf, depth, height, width, min_intensity) = min_intensity *sum([psf(x -width/2, y - height/2, z - depth/2) for x in 1:width, y in 1:height, z in 1:depth] .^2)


"""
    get_conv_cand_dot_coords(image, psf, thresh)

Convolve the image and PSF, then return the coordinates of the result above the threshold and candidate
dot coordinates as dataframe
"""
function get_conv_cand_dot_coords(imagestack, psf, psf_cutoff, min_peak)
    println("size(imagestack): ", size(imagestack))
    #psf_grid = [psf(i,j,k) for i in -20:20, j in -20:20, k in -20:20]# .> psf_cutoff
    psf_grid = [psf(i,j) for i in -20:20, j in -20:20]# .> psf_cutoff

    psf_grid[psf_grid .< psf_cutoff] .= 0
    #thresh = sum(psf_grid .^2) * min_intensity
    #zproj = reshape(sum(sum(psf_grid, dims=1), dims=2), 41)
    #yproj = reshape(sum(sum(psf_grid, dims=1), dims=3), 41)
    #xproj = reshape(sum(sum(psf_grid, dims=3), dims=1), 41)
    col_proj = reshape(sum(psf_grid, dims=1), 41)
    row_proj = reshape(sum(psf_grid, dims=2), 41)


    psf_ncols = findlast(v -> v > 0, col_proj) - findfirst(v -> v > 0, col_proj) + 1
    psf_nrows = findlast(v -> v > 0, row_proj) - findfirst(v -> v > 0, row_proj) + 1
    #psf_depth = findlast(v -> v > 0, zproj) - findfirst(v -> v > 0, zproj) + 1

    #thresh =  get_conv_threshold(psf, psf_depth, psf_height, psf_width, min_peak)
    thresh =  get_conv_threshold(psf, 1, psf_nrows, psf_ncols, min_peak)



    #im_depth, im_height, im_width = size(imagestack)
    
    im_nrows, im_ncols = size(imagestack)

    cand_dot_coords = []
    #observed = zeros(im_depth*im_height*im_width)
    #psf_mod = zeros(im_depth*im_height*im_width)
    println("begin convolutions")
    #for z1 in 1:0.5:im_depth, y1 in 1:0.5:im_height, x1 in 1:0.5:im_width
    #for z1 in 1:1:im_depth, y1 in 1:1:im_height, x1 in 1:1:im_width
    for r1 in 1:1:im_nrows, c1 in 1:1:im_ncols
        conv_val = 0
        #observed .= 0.0
        #psf_mod .= 0.0
        #i = 1
        #for z2 in maximum([1, floor(Int64, z1 - psf_depth/2)]):minimum([im_depth, ceil(Int64, z1 + psf_depth/2)])
        for r2 in maximum([1, floor(Int64, r1 - psf_nrows/2)]):minimum([im_nrows, ceil(Int64, r1 + psf_nrows/2)])
            for c2 in maximum([1, floor(Int64, c1 - psf_ncols/2)]):minimum([im_ncols, ceil(Int64, c1 + psf_ncols/2)])
                #conv_val += imagestack[z2, y2, x2] * psf(x1-x2, y2-y2, z1-z2)
                conv_val += imagestack[r2, c2] * psf(r1-r2, c1-c2)
                #observed[i] = imagestack[z2, y2, x2]
                #psf_mod[i] = psf(x1-x2, y2-y2, z1-z2)
            end
        end
        #end
        if conv_val > thresh
            #push!(cand_dot_coords, [z1, y1, x1])
            push!(cand_dot_coords, (r1, c1))
        end
    end
    return cand_dot_coords #, (r1, c1)
end

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

sparse_psf = DataFrame(Dict("row" => psf_is, "col" => psf_js, "v" => psf_vs))



function get_sparse_dot_vec(dot_row, dot_col, row_dim, col_dim, _sparse_psf_df, j)
    sparse_dot_rep = copy(_sparse_psf_df)
            
    # move to dot location
    sparse_dot_rep[!, "row"] .+= dot_row
    sparse_dot_rep[!, "col"] .+= dot_col
    
    # remove out of bound pixels
    filter!(p -> (p.row > 0) && (p.row <= row_dim) && (p.col > 0) && (p.col <= col_dim), sparse_dot_rep)

    # convert row-col of sparse block-pseudocolor image represntation to i of sparse design matrix
    is = sparse_dot_rep.row .+ row_dim*(sparse_dot_rep.col .- 1)
    return DataFrame(Dict("row" => is, "col" => fill(j, length(is)), "v" => sparse_dot_rep.v))
end

function reindex_sprep(img, A_csc_df)

    ys = reshape(img, prod(size(img)))
    sorted_unique_rows = sort(unique(A_csc_df.row))
    ydf = DataFrame(Dict("row" => sorted_unique_rows, "y" => ys[sorted_unique_rows]))
        
    reindex_row_dict = Dict([(row0 => rowp) for (rowp, row0) in enumerate(sorted_unique_rows)])
    reindex_row(row) = reindex_row_dict[row]
    A_csc_df[!,"row"] .= reindex_row.(A_csc_df[!,"row"])
    ydf[!,"row"] .= reindex_row.(ydf[!,"row"])

    filter!(r -> r.y > 0, ydf)
    ydf = sort(unique(ydf))
    return ydf, A_csc_df
end


#for i in 1:nbcrounds
   # hyperstack[i, 1, :, :] = imstack[:,:,i]
#end

img = reinterpret.(UInt16, channelview(load(snakemake.input[1])))
img = round.(UInt16, img)


#img = hyperstack[:, :, s]

get_sparse_dot_vec(dot_row, dot_col, j) = get_sparse_dot_vec(dot_row, dot_col, size(img)[1], size(img)[2], sparse_psf, j)
get_sparse_dot_vec(t :: Tuple{Int64, Tuple{Int64, Int64}}) = get_sparse_dot_vec(t[2][1], t[2][2], t[1])


#hyperstack = reshape(hyperstack, vcat([1], collect(size(hyperstack)))...)


ccoords = get_conv_cand_dot_coords(img, psf, psf_cutoff, min_peak)

A_csc_df = vcat(get_sparse_dot_vec.(enumerate(ccoords))...)

ydf, A_csc_df_reindex = reindex_sprep(img, A_csc_df)

A_spn = sparse(UInt32.(A_csc_df.row), UInt32.(A_csc_df.col), Float64.(A_csc_df.v))

constraints = zeros(2, maximum(A_csc_df.col))
constraints[2,:] .= 2^16

Y = SparseVector(maximum(A_csc_df.row), UInt32.(ydf.row), UInt16.(ydf.y))

lasso_path = glmnet(A_spn, Y, constraints = constraints)

#choose the lambda value that produces the smallest non-negative intercept
intercepts = lasso_path.a0
intercepts[intercepts .< 0] .= Inf
lambda_ind = argmin(intercepts)
selected_dot_inds = lasso_path.betas[:, lambda_ind] .> min_peak
selected_ccoords = ccoords[selected_dot_inds]
cand_dots = DataFrame(hcat(collect.(selected_ccoords)...)', ["x", "y"])

cand_dots[!, "z"] .= 1
cand_dots[!, "s"] .= 0
cand_dots[!, "s_xy"] .= 0
cand_dots[!, "s_z"] .= 0
cand_dots[!, "w"] .= 1
cand_dots[!, "block"] .= r
cand_dots[!, "pseudocolor"] .= pc
cand_dots[!, "hyb"] .= hyb
write_parquet(snakemake.output[1], cand_dots)