using DataFrames, CSV, Images, FileIO

lbld_img = reinterpret.(UInt8, channelview(load(snakemake.input[1])))

cell_nums = sort(filter(cn -> cn != 0, unique(lbld_img)))
cell_xmaxes = []
cell_xmins = []
cell_ymaxes = []
cell_ymins = []
for cell_num in cell_nums
    push!(cell_xmaxes, maximum(map(c -> c[1], findall(s -> s == 1, maximum(lbld_img .== cell_num, dims=2)))))
    push!(cell_xmins, minimum(map(c -> c[1], findall(s -> s == 1, maximum(lbld_img .== cell_num, dims=2)))))
    push!(cell_ymaxes, maximum(map(c -> c[2], findall(s -> s == 1, maximum(lbld_img .== cell_num, dims=1)))))
    push!(cell_ymins, minimum(map(c -> c[2], findall(s -> s == 1, maximum(lbld_img .== cell_num, dims=1)))))
end

df = DataFrame(:cell_num => cell_nums, :xmin => cell_xmins, :xmax => cell_xmaxes, :ymin => cell_ymins, :ymax => cell_ymaxes)

CSV.write(snakemake.output[1], df)