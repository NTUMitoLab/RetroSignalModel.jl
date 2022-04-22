using RetroSignalModel
using ModelingToolkit
using DataFrames
using Optim
using CSV

optimoptions = Optim.Options(iterations=10^5, show_trace=true, show_every=1000)

parammaps = map(1:10) do iter
    @show iter
    _, parammap = optim_params(; optimoptions)
    parammap
end

parammapssymbol = map(parammaps) do x
    Dict(Symbol(k) => v for (k, v) in x)
end

df = DataFrame(parammapssymbol)

CSV.write("solution_rtgMTK_optim.csv", df)
