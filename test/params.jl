# Test scan_params
const nruns = 10

sim = scan_params(; trajectories=10, saveall=true)
@test length(sim) == nruns

# Test optim_params
import Optim: Options, iterations
res = optim_params(; optimoptions=Options(iterations=nruns, show_trace=true))

@test iterations(res) == nruns
