# Parameter searching to meet the conditions
using SciMLBase
using SteadyStateDiffEq
using OrdinaryDiffEq
using ModelingToolkit
using Optim

# Distributions of RTG1 and RTG3 proteins
rtg13_nucleus(sol) = sol[Rtg13I_n] + sol[Rtg13A_n]
rtg13_cytosol(sol) = sol[Rtg13I_c] + sol[Rtg13A_c]
rtg3_nucleus(sol) = rtg13_nucleus(sol) + sol[Rtg3A_n] + sol[Rtg3I_n]
rtg3_cytosol(sol) = rtg13_cytosol(sol) + sol[Rtg3A_c] + sol[Rtg3I_c]
rtg1_nucleus(sol) = rtg13_nucleus(sol) + sol[Rtg1_n]
rtg1_cytosol(sol) = rtg13_cytosol(sol) + sol[Rtg1_c]

"""
Randomly scan for parameters that satisfy the boolean conditions in the retrograde (RTG) signalling model.
"""
function scan_params(
    Model=RtgMTK;
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    knockoutlevel=1e-6,
    nuclearRatioThreshold=10,
    trajectories=1000,
    batch_size=trajectories,
    proteinlevels=STRESSED,
    steadyStateSolver=DynamicSS(Rodas5()),
    ensembleSolver=EnsembleThreads(),
    ntarget=100,
    saveall=false,
    rollparams=() -> exp10(6 * (rand() - 0.5)),
    rollhill=() -> rand(1.0:0.5:5.0))
    # Boolean conditions for nuclear accumulation
    conds = load_conditions(datafile)

    @named sys = Model(ONE_SIGNAL; proteinlevels, recordNruns=true)

    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
    idxnRuns = param2idx[nRuns]
    idxnS = param2idx[n_S]

    prob = SteadyStateProblem(sys, resting_u0(sys))

    # Initalize parameters for the first batch
    batchParams = zeros(batch_size, length(parameters(sys)))

    # Populate parameters for the batches
    function populate_batch_params!(batchParams)
        for i in eachindex(batchParams)
            batchParams[i] = rollparams()
        end

        for i in 1:batch_size
            batchParams[i, idxnS] = rollhill()
        end
    end

    populate_batch_params!(batchParams)

    # select random param for a fresh problem
    # repeat for each conditions
    function prob_func(prob, i, iter)
        params = prob.p
        params .= view(batchParams, i % batch_size + 1, :)
        params[idxnRuns] = iter

        # Adjust params according to conditions
        cond = conds[iter]
        params[param2idx[ΣRtg1]] = ifelse(cond[:Rtg1] == 0, knockoutlevel, proteinlevels[ΣRtg1])
        params[param2idx[ΣRtg2]] = ifelse(cond[:Rtg2] == 0, knockoutlevel, proteinlevels[ΣRtg2])
        params[param2idx[ΣRtg3]] = ifelse(cond[:Rtg3] == 0, knockoutlevel, proteinlevels[ΣRtg3])
        params[param2idx[ΣMks]] = ifelse(cond[:Mks] == 0, knockoutlevel, proteinlevels[ΣMks])
        params[param2idx[mul_S]] = cond[:s]

        remake(prob, p=params)
    end

    function pass_cond(idx, conds, sol)
        cond = conds[idx]
        passed = false
        if cond[:gfp] == "rtg3"
            if cond[:Trans2Nuc] == 1
                passed = rtg3_nucleus(sol) > nuclearRatioThreshold * rtg3_cytosol(sol)
            else
                passed = nuclearRatioThreshold * rtg3_nucleus(sol) < rtg3_cytosol(sol)
            end
        elseif cond[:gfp] == "rtg1"
            if cond[:Trans2Nuc] == 1
                passed = rtg1_nucleus(sol) > nuclearRatioThreshold * rtg1_cytosol(sol)
            else
                passed = nuclearRatioThreshold * rtg1_nucleus(sol) < rtg1_cytosol(sol)
            end
        end
        return passed
    end

    # rerun when passed a condition and number for runs < number of conditions
    # otherwise, do not rerun
    function output_func(sol, i)
        idx = Int(sol.prob.p[idxnRuns])
        rerun = (idx < length(conds)) && pass_cond(idx, conds, sol)
        return (sol, rerun)
    end

    # Only add tests that passed all conditions
    function reduction(u, batch, I)

        for sol in batch
            idx = Int(sol.prob.p[idxnRuns])
            if saveall || (idx == length(conds) && pass_cond(idx, conds, sol))
                u = push!(u, sol)
            end
        end

        populate_batch_params!(batchParams)

        return (u, length(u) >= ntarget)
    end

    ensprob = EnsembleProblem(prob; output_func, prob_func, reduction)

    sim = solve(ensprob, steadyStateSolver, ensembleSolver; trajectories, batch_size)
end

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
    optimoptions=Optim.Options(iterations=10^6))

    # Boolean conditions
    conds = load_conditions(datafile)

    @named sys = Model(ONE_SIGNAL; proteinlevels)
    prob = SteadyStateProblem(sys, resting_u0(sys))

    # Mapping between model parameters and the cost function inputs
    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
    idxnS = param2idx[n_S]
    idxΣRtg1 = param2idx[ΣRtg1]
    idxΣRtg2 = param2idx[ΣRtg2]
    idxΣRtg3 = param2idx[ΣRtg3]
    idxΣMks = param2idx[ΣMks]
    idxmul_S = param2idx[mul_S]

    xidx2params = [i for (k, i) in param2idx if !any(isequal(k), (ΣRtg1, ΣRtg2, ΣRtg3, ΣMks, mul_S))]
    xidxnS = findfirst(isequal(idxnS), xidx2params)

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

            for i in 1:length(x)
                params[xidx2params[i]] = x[i]
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
            count += score
        end

        return count / length(conds)
    end

    x0 = similar(xidx2params, Float64)
    x0 .= xinit
    lb = similar(x0)
    lb .= lowerbound
    lb[xidxnS] = hilllowerbound
    ub = similar(x0)
    ub .= upperbound
    ub[xidxnS] = hillupperbound

    res = Optim.optimize(cost, lb, ub, x0, optimsolver, optimoptions)
end
