using SeqFISH_ADCG
using Images
using FileIO
using DataFrames
using CSV

fov_width = snakemake.params["fov_width"]
sigma = snakemake.params["sigma"]
min_weight = snakemake.params["min_weight"] #200.0
max_iters = 400
max_cd_iters = 40
final_loss_improvement = 500.0

gblur = GaussBlur2D(sigma, sigma, fov_width)

img = reinterpret.(UInt16, channelview(load(snakemake.input[1])))

inputs = (img, sigma, sigma, 0.0, 0.0, final_loss_improvement, min_weight, max_iters, max_cd_iters, "ADCG")
record = SeqFISH_ADCG.fit_tile(inputs)

pnts_df = record.last_iteration

#if nrow(points) == 0
#    pnts_df = DataFrame(x=[],y=[],s=[],w=[])
#else
#    pnts_df = DataFrame(points', [:x,:y, :s, :w])
#end
pnts_df[!, "cellid"] .= 1

CSV.write(snakemake.output[1], pnts_df)

