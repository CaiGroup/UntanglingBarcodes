using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames
using CSV
using Gurobi
using JuMP

cb = DataFrame(CSV.File("resources/q11_n10_k7_w4cb_assign_3t3_ribo_542_prethinned_attempt1.csv"))

cb_no_genes = select(cb, Not(:gene))

nblocks = ncol(cb_no_genes)# - 1
npseudocolors = maximum(Matrix(cb_no_genes))

println("npseudocolors: $npseudocolors")
println("nblocks: $nblocks")

binary_cb = BitArray(zeros(nrow(cb_no_genes), npseudocolors*nblocks))

cb_channel_table = DataFrame(zeros(Int64, npseudocolors*nblocks, 3), [:block, :pseudocolor, :channel])

for (i, cw) in enumerate(eachrow(cb_no_genes))
    for (r, pc) in enumerate(cw)
        if pc != 0
            #println("pc: ", pc)
            #println(npseudocolors*(r-1))
            binary_cb[i, pc + npseudocolors*(r-1)] = 1
            cb_channel_table[pc + npseudocolors*(r-1),"block"] = r
            cb_channel_table[pc + npseudocolors*(r-1),"pseudocolor"] = pc
        end
    end
end

model = Model(Gurobi.Optimizer)

MOI.set(model, MOI.TimeLimitSec(), 3600*2)

@variable(model, channel_indicator[1:(nblocks*npseudocolors), 1:3], Bin)

@constraint(model, sum(channel_indicator, dims=2) .== 1)
@constraint(model, sum(channel_indicator, dims=1) .<= ceil(nblocks*npseudocolors/3))


cw_ch_ro_cnt = binary_cb * channel_indicator

@variable(model, ro_ch_cnt_ind[1:(nrow(cb)), 1:3, 1:4], Bin)

@constraint(model, sum(ro_ch_cnt_ind, dims=3) .<= 1)
@constraint(model, cw_ch_ro_cnt .== sum(ro_ch_cnt_ind .* reshape([1,2,3,4],1,1,4), dims=3))

#@constraint(model, ro_ch_max .>= cw_ch_ro_cnt)
#@constraint(model, .- ro_ch_diff .<=  .- cw_ch_ro_cnt)
#@objective(model, Max, sum(cw_ch_ro_cnt .^ 2))
@objective(model, Max, sum(ro_ch_cnt_ind[:,:,3:4]))
optimize!(model)

ch_ind_res = round.(Int64, value.(channel_indicator))

cb_channel_table[!, "channel"] .= maximum((ch_ind_res' .* [1,2,3])', dims = 2)

CSV.write("block_pseudocolor_channel_time2hr_linear.csv", cb_channel_table)