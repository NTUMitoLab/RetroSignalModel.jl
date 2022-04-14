using SciMLBase

# Test find_steady_states
res = find_steady_states(trajectories=10)

@test res isa EnsembleSolution
@test length(res) > 0
