using RetroSignalModel
using ModelingToolkit
using DataFrames
using Optim
using CSV

versioninfo()

optimoptions = Optim.Options(iterations=10^5)

res = optim_params_threads(120)

parammaps = map(res) do x
    Dict(Symbol(k) => v for (k, v) in x.parammap)
end

CSV.write("solution_rtgMTK_optim.csv", DataFrame(parammaps))
