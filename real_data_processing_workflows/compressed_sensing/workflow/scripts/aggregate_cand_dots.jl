using DataFrames, SeqFISHSyndromeDecoding, Parquet

#df = vcat(DataFrame.(CSV.File(values(snakemake.input)))...)
df = vcat(DataFrame.(read_parquet.([snakemake.input[i] for i in 1: length(snakemake.input)]))...)


df[!, "pseudocolor"] .%= 20
#df[!, "z"] = zeros(nrow(df))

SeqFISHSyndromeDecoding.sort_readouts!(df)

write_parquet(snakemake.output[1], df)