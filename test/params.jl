# Test scan_params
const nruns = 10
nthreads = Threads.nthreads()

# Test optim_params
using Optim
res, _ = optim_params(; optimoptions=Optim.Options(iterations=nruns, show_trace=true))

@test Optim.iterations(res) == nruns

ensres = optim_params_threads(nthreads, optimoptions=Optim.Options(iterations=nruns))

@test length(ensres) == nthreads
