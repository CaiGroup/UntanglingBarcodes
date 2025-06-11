using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames
using CSV
using JuMP
using Gurobi
#addprocs(20)
cb = DataFrame(CSV.File("resources/q11_n10_k7_w4cb_assign_3t3_ribo_542_prethinned_attempt1.csv"))
sort!(cb, :gene)
genes = DataFrame(CSV.File("resources/ribo&HE_genes_100kfpkm.csv"))
sort!(genes, :tracking_id)


cb_no_genes = select(cb, Not(:gene))

println("nrow(cb): ", nrow(cb))

genes = DataFrame(CSV.File("resources/ribo&HE_genes_100kfpkm.csv"))
println("nrow(fpkm): ", nrow(genes))


npseudocolors = maximum(Array(cb_no_genes))
nblocks = ncol(cb)# - 1
println("npseudocolors: $npseudocolors")
println("nblocks: $nblocks")

binary_cb = BitArray(zeros(nrow(cb_no_genes), npseudocolors*nblocks))

for (i, cw) in enumerate(eachrow(cb_no_genes))
    for (r, pc) in enumerate(cw)
        if pc != 0
            #println("pc: ", pc)
            #println(npseudocolors*(r-1))
            binary_cb[i, pc + npseudocolors*(r-1)] = 1
        end
    end
end

#@assert length(unique(sum(binary_cb, dims=1))) == 1

model = Model(Gurobi.Optimizer)

MOI.set(model, MOI.TimeLimitSec(), 40)##70000) # set time limit


@variable(model, pool1_indicator[1:nrow(genes)], Bin)
println("add pool1")
@constraint(model, sum(pool1_indicator) == nrow(genes)/2)


println("get hyb_fpkm")
hyb_fpkm_pool1 = (genes[:, "3T3 B1"] .* pool1_indicator)'*binary_cb
hyb_fpkm_pool2 = (genes[:, "3T3 B1"] .* (1 .- pool1_indicator))'*binary_cb


@variable(model, mx)

println("adding max mn constraints")

@constraint(model, hyb_fpkm_pool1 .<= mx)
@constraint(model, hyb_fpkm_pool2 .<= mx)


println("optimizing")
@objective(model, Min, mx)
optimize!(model)

println(value(mx))

println("done optimizing")

#assign_order = (1:nrow(genes))'*value.(gene_cw_assign_indicator)

int_pool1_indicator = Int64.(value.(pool1_indicator))
pool1 = cb[findall(x -> x == 1, int_pool1_indicator), :]
pool2 = cb[findall(x -> x == 0, int_pool1_indicator), :]


CSV.write("q11_n10_k7_w4cb_assign_3t3_ribo_542_prethinned_pool1_attempt1.csv", pool1)
CSV.write("q11_n10_k7_w4cb_assign_3t3_ribo_542_prethinned_pool2_attempt1.csv", pool2)