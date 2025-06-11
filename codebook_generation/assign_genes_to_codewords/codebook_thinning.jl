using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DataFrames
using CSV
using JuMP
using Gurobi

cb = DataFrame(CSV.File("resources/RS_q11_k7_w4cb.csv"))

ncws = 542

#get distance matrix
D = zeros(Int64, nrow(cb), nrow(cb))



model = Model(Gurobi.Optimizer)

MOI.set(model, MOI.TimeLimitSec(), 3600 * 16) # set time limit to 12 hours


@variable(model, X[1:nrow(cb), 1:nrow(cb)], Bin, Symmetric)

@variable(model, chosen_cws[1:nrow(cb)], Bin)

for (i, ri) in enumerate(eachrow(cb)), (j, rj) in enumerate(eachrow(cb))
    D[i, j] = sum(Array(ri) .!= Array(rj))
    @constraint(model, X[i,j] <= chosen_cws[i])
    @constraint(model, X[i,j] <= chosen_cws[j])
    @constraint(model, X[i,j] >= chosen_cws[i] + chosen_cws[j]-1)
end


nrounds = ncol(cb)# - 1
npseudocolors = maximum(Matrix(cb))

println("npseudocolors: $npseudocolors")
println("nrounds: $nrounds")

binary_cb = BitArray(zeros(nrow(cb), npseudocolors*nrounds))


for (i, cw) in enumerate(eachrow(cb))
    for (r, pc) in enumerate(cw)
        if pc != 0
            #println("pc: ", pc)
            #println(npseudocolors*(r-1))
            binary_cb[i, pc + npseudocolors*(r-1)] = 1
        end
    end
end


@constraint(model, sum(chosen_cws) == ncws)
@constraint(model, chosen_cws[1] == 1)

cws_per_hyb = chosen_cws' * binary_cb

@constraint(model, cws_per_hyb .>= ncws*4/100-5)
@constraint(model, cws_per_hyb .<= ncws*4/100+5)

@objective(model, Max, sum(X .* D))

optimize!(model)

thinned_cb = cb[round.(Bool,value.(chosen_cws)),:]

CSV.write("thinned_cb_542cws_even_per_hyb_16hr.csv", thinned_cb)
