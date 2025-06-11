using Pkg
Pkg.activate(".")
Pkg.instantiate()


using DataFrames
using CSV
using JuMP
using Gurobi
using Distributed
#addprocs(20)
cb = DataFrame(CSV.File("resources/thinned_cb_542cws_try2.csv"))

println("nrow(cb): ", nrow(cb))

genes = DataFrame(CSV.File("resources/ribo&HE_genes_100kfpkm.csv"))
println("nrow(fpkm): ", nrow(genes))
sort!(genes, "3T3 B1", rev=true)


npseudocolors = maximum(Array(cb))
nblocks = ncol(cb)# - 1
println("npseudocolors: $npseudocolors")
println("nblocks: $nblocks")

binary_cb = BitArray(zeros(nrow(cb), npseudocolors*nblocks))

for (i, cw) in enumerate(eachrow(cb))
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

MOI.set(model, MOI.TimeLimitSec(), 4000)##70000) # set time limit to 12 hours

@variable(model, gene_cw_assign_indicator[1:nrow(genes), 1:nrow(cb)], Bin)
println("add gene_cw_assign_indicator constraint 1")
@constraint(model, sum(gene_cw_assign_indicator, dims=1) .<= 1)
println("add gene_cw_assign_indicator constraint 2")
@constraint(model, sum(gene_cw_assign_indicator, dims=2) .== 1)

println("get hyb_fpkm")
hyb_fpkm = genes[:, "3T3 B1"]'*gene_cw_assign_indicator*binary_cb


@variable(model, mx)

println("adding max mn constraints")

@constraint(model, hyb_fpkm .<= mx)

println("optimizing")
@objective(model, Min, mx)
optimize!(model)

println(value(mx))

println("done optimizing")

assign_order = (1:nrow(genes))'*value.(gene_cw_assign_indicator)

insertcols!(cb, 1, :gene_num => reshape(Int64.(assign_order), nrow(cb)))

filter!(r -> r.gene_num != 0, cb)
sort!(cb, "gene_num")

cb = select(cb, Not(:gene_num))

insertcols!(cb, 1, :gene => genes.tracking_id)

CSV.write("q11_n10_k7_w4cb_assign_3t3_ribo_542_prethinned_attempt1.csv", cb)