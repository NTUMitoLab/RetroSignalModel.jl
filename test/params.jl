using Optim
using DataFrames

const nruns = 10
nthreads = Threads.nthreads()

res, params = optim_params(; optimoptions=Optim.Options(iterations=nruns, show_trace=false))
@test Optim.iterations(res) == nruns

df = optim_params_threads(nthreads, optimoptions=Optim.Options(iterations=nruns))
@test df isa DataFrame
