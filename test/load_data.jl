import RetroSignalModel as rs

conds = load_conditions()

const condkeys = (:Rtg1, :Rtg2, :Rtg3, :Mks, :s, :gfp, :Trans2Nuc)

@testset "Name availability (Conditions): $k" for cond in conds, k in condkeys
    @test k in keys(cond)
end

const paramkeys = (
    rs.n_S,
    rs.ksV,
    rs.ksD,
    rs.k2I,
    rs.k2M,
    rs.kn2M,
    rs.kBM,
    rs.knBM,
    rs.k13I,
    rs.k13IV,
    rs.k13ID,
    rs.k3A_c,
    rs.k3I_c,
    rs.k3I_n,
    rs.k13_c,
    rs.kn13_c,
    rs.k13_n,
    rs.kn13_n,
    rs.k1in,
    rs.k1out,
    rs.k3inA,
    rs.k3outA,
    rs.k3inI,
    rs.k3outI)

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
