using SciMLBase
using SteadyStateDiffEq

# Test find_steady_states
res = find_steady_states(trajectories=10)

@test res isa EnsembleSolution
@test res[1] isa SteadyStateSolution
@test length(res) > 0
