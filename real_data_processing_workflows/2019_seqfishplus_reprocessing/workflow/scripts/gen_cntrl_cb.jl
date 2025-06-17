using DelimitedFiles
using CSV
using DataFrames
using LinearAlgebra
using Random


gene_only_cb_fname = snakemake.input[1]
cb_df = DataFrame(CSV.File(gene_only_cb_fname))
cb_df.gene = String.(cb_df.gene)

cb = Matrix(UInt8.(cb_df[!,2:end]))
q = maximum(cb)
cb = cb .% q

H = readdlm(snakemake.input[2], Int64)

n_genes_encoded, ncols = size(cb)


cws = []

cws = [cb[i,:] for i in 1:n_genes_encoded]
sort!(cws)

control_cws = []
for i = 0:(q-1), j = 0:(q-1), k = 0:(q-1), l = 0:(q-1)
    cw = [i, j, k, l]
    if [i j k l] ⋅ H % q == 0
        loc = searchsorted(cws, cw)
            #println("cw: ", cw)
            #if cw ∈ cws && isempty(loc)
                #println(true)
            #end
        if isempty(loc)
            push!(control_cws, cw)
        end
    end
end

if all(H .== 0)
    @assert length(cws) + length(control_cws) == Int64(q)^4
else
    @assert length(cws) + length(control_cws) == Int64(q)^3
end

# if not using all control codewords, take a random subset
if (snakemake.params["n_neg_cntrl_cws"] != "all") && (snakemake.params["n_neg_cntrl_cws"] < length(control_cws))
    # draw codewords without replacement
    control_cws = sort(shuffle(control_cws)[1:snakemake.params["n_neg_cntrl_cws"]])
end

ctrl_cb = DataFrame(zeros(UInt8, length(control_cws), 4), ["block1", "block2", "block3", "block4"])


for i = 1:nrow(ctrl_cb)
    ctrl_cb[i,:] = pop!(control_cws)
end

insertcols!(ctrl_cb,1,"gene" => fill("negative_control", nrow(ctrl_cb)))
append!(cb_df, ctrl_cb)
#rename!(cb_df,"gene name"=>"gene")
#rename!(cb_df,:hyb1=>:block1,:hyb2=>:block2,:hyb3=>:block3,:hyb4=>:block4)

cb_df[!, 2:end] .= cb_df[!, 2:end] .% q

CSV.write(snakemake.output[1], cb_df)
