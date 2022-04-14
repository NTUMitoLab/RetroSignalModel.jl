# Test scan_params
const nruns = 10

# Test optim_params
using Optim
res, _ = optim_params(; optimoptions=Optim.Options(iterations=nruns, show_trace=true))

@test Optim.iterations(res) == nruns
