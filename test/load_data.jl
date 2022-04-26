conds = load_conditions()

const condkeys = (:Rtg1, :Rtg2, :Rtg3, :Mks, :s, :gfp, :Trans2Nuc)

@testset "Name availability (Conditions): $k" for cond in conds, k in condkeys
    @test k in keys(cond)
end

const paramkeys = (
    :n_S,
    :ksV,
    :ksD,
    :k2I,
    :k2M,
    :kn2M,
    :kBM,
    :knBM,
    :k13I,
    :k13IV,
    :k13ID,
    :k3A_c,
    :k3I_c,
    :k3I_n,
    :k13_c,
    :kn13_c,
    :k13_n,
    :kn13_n,
    :k1in,
    :k1out,
    :k3inA,
    :k3outA,
    :k3inI,
    :k3outI)

optimparams = load_parameters("solution_rtgMTK_optim.csv")

# Solutions of rtgM4
@testset "Name availability (Optim Parameters): $k" for param in optimparams, k in paramkeys
    @test k in keys(param)
end

randparams = load_parameters("solution_rtgM4.csv")

# Solutions of rtgM4
@testset "Name availability (Random Parameters): $k" for param in randparams, k in paramkeys
    @test k in keys(param)
end
