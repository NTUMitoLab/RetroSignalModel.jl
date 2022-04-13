using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq

import ModelingToolkit as mtk

"""
Find steady states of the ODE system `sys` by randomly
"""
function find_steady_states(sys;)

    statesmap = Dict(k => i for (i, k) in enumerate(mtk.states(sys)))

    # Create a set of initial conditions respecting conservation relationships for the ensemble
    function make_u0!(u0)

        # Randomly split total into `n` fractions
        function _frac(n::Int, total=1.0)
            total .* diff([zero(total); sort([rand(typeof(total)) for _ in 1:(n-1)]); one(total)])
        end

        # Split ΣMks into 3
        frac = _frac(3, ΣMks)
        u0[statesmap[Mks]] = frac[1]
        u0[statesmap[BmhMks]] = frac[2]
        u0[statesmap[Rtg2Mks_c]] = frac[3]
        u0[statesmap[Bmh]] = Σbmh - u0[statesmap[BmhMks]]

        frac = _frac(2, ΣRtg2 - u0[statesmap[Rtg2Mks_c]])
        u0[statesmap[Rtg2_ina_c]] = frac[1]
        u0[statesmap[Rtg2_act_c]] = frac[2]

        # Split ΣRtg1 into 6
        frac = _frac(6, ΣRtg1)
        u0[statesmap[Rtg13_a_c]] = frac[1]
        u0[statesmap[Rtg13_i_c]] = frac[2]
        u0[statesmap[Rtg1_c]] = frac[3]
        u0[statesmap[Rtg1_n]] = frac[4]
        u0[statesmap[Rtg13_a_n]] = frac[5]
        u0[statesmap[Rtg13_i_n]] = frac[6]

        # Split ΣRtg3 not in Rtg13 complex into 4
        frac = _frac(4, ΣRtg3 - u0[statesmap[Rtg13A_c]] - u0[statesmap[Rtg13I_c]] - u0[statesmap[Rtg13A_n]] - u0[statesmap[Rtg13I_n]])

        u0[statesmap[Rtg3_i_c]] = frac[1]
        u0[statesmap[Rtg3_a_c]] = frac[2]
        u0[statesmap[Rtg3_a_n]] = frac[3]
        u0[statesmap[Rtg3_i_n]] = frac[4]

        return u0
    end

    function prob_func(prob, i, repeat)
        u0 = make_u0!(prob.u0, Σbmh, ΣMks, ΣRtg1, ΣRtg2, ΣRtg3)
        remake(prob, u0=u0)
    end



end


function prob_func(prob, i, repeat)
    make_u0!(prob.u0, Σbmh, ΣMks, ΣRtg1, ΣRtg2, ΣRtg3)
    prob
end

"""Reject invalid (NaNs and negatives) and duplicate results"""
function reduction(u, batch, I; negtol=-1e-6)
    for result in batch
        # Skip invalid (NaNs and negatives) results
        if any(isnan.(result)) || any(result .< negtol)
            continue
        end

        result .= max.(0, result)
        isUnique = true

        # Only save unique solutions
        for existing in u
            if all(isapprox.(existing, result, rtol=1e-4))
                isUnique = false
                break
            end
        end

        if isUnique
            u = push!(u, result)
        end
    end
    return (u, false)
end
