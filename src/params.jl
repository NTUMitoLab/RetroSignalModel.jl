# Parameter searching to meet the conditions
using SciMLBase
using SteadyStateDiffEq
using OrdinaryDiffEq
using ModelingToolkit
using Optim

"""
Find parameters that satisfy the boolean conditions in the retrograde (RTG) signalling model using Optimization methods.
"""
function optim_params(
    Model=RtgMTK;
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    knockoutlevel=1e-4,
    proteinlevels=STRESSED,
    steadyStateSolver=DynamicSS(Rodas5()),
    lowerbound=1e-3,
    upperbound=1e3,
    hilllowerbound=1.0,
    hillupperbound=7.0,
    xinit=1.0,
    optimsolver=Optim.SAMIN(),
    optimoptions=Optim.Options(iterations=10^5, show_trace=true, show_every=1000),
    targetratio=10)

    # Boolean conditions
    conds = load_conditions(datafile)

    @named sys = Model(ONE_SIGNAL; proteinlevels)
    prob = SteadyStateProblem(sys, resting_u0(sys))

    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
    idxΣRtg1 = param2idx[ΣRtg1]
    idxΣRtg2 = param2idx[ΣRtg2]
    idxΣRtg3 = param2idx[ΣRtg3]
    idxΣMks = param2idx[ΣMks]
    idxmul_S = param2idx[mul_S]

    params_optim = [k for k in parameters(sys) if !any(isequal(k), (ΣRtg1, ΣRtg2, ΣRtg3, ΣMks, ΣBmh, mul_S))]

    # Mapping indices of x in Optim to param indices in the ODE system
    xi2params = [param2idx[k] for k in params_optim]
    idxnS = param2idx[n_S]
    xidxnS = findfirst(isequal(idxnS), xi2params)

    # Initial conditions and lower / upper bounds for Optim
    x0 = similar(xi2params, Float64)
    x0 .= xinit
    lb = similar(x0)
    lb .= lowerbound
    lb[xidxnS] = hilllowerbound
    ub = similar(x0)
    ub .= upperbound
    ub[xidxnS] = hillupperbound

    # The cost is bounded by target protein concentration ratio
    scorecap = -log10(targetratio)

    # Cost function
    function cost(x)
        count = 0.0
        for cond in conds
            params = copy(prob.p)

            # Adjust params according to conditions
            params[idxΣRtg1] = ifelse(cond[:Rtg1] == 0, knockoutlevel, proteinlevels[ΣRtg1])
            params[idxΣRtg2] = ifelse(cond[:Rtg2] == 0, knockoutlevel, proteinlevels[ΣRtg2])
            params[idxΣRtg3] = ifelse(cond[:Rtg3] == 0, knockoutlevel, proteinlevels[ΣRtg3])
            params[idxΣMks] = ifelse(cond[:Mks] == 0, knockoutlevel, proteinlevels[ΣMks])
            params[idxmul_S] = cond[:s]

            # Align parameters to the vector to be optimized
            for i in 1:length(x)
                params[xi2params[i]] = x[i]
            end

            # Calculate cost for this steady state solution
            # Since the objective is nuclear accumulation or not
            # The cost is the logarithm of protein concentration ratios (nuclear vs cytosol)
            sol = solve(remake(prob, p=params), steadyStateSolver)

            if cond[:gfp] == "rtg3"
                score = -log10(rtg3_nucleus(sol) / rtg3_cytosol(sol)) * ifelse(cond[:Trans2Nuc] == 1, 1, -1)
            elseif cond[:gfp] == "rtg1"
                score = -log10(rtg1_nucleus(sol) / rtg1_cytosol(sol)) * ifelse(cond[:Trans2Nuc] == 1, 1, -1)
            else
                score = 0.0
            end
            score = max(score, scorecap)
            count += score
        end

        return count / length(conds)
    end
    res = Optim.optimize(cost, lb, ub, x0, optimsolver, optimoptions)
    parammap = Dict(params_optim .=> Optim.minimizer(res))

    return (res=res, parammap=parammap)
end
