
using DrWatson
@quickactivate "Floating_Point_Error_Reproducer"

using DifferentialEquations
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using LineSearches

include("Structs.jl")
include("EOMs.jl")
include("CBs.jl")
include("CR3BPIMF_ShootingFunction_STM.jl")

function main()

    # Load CoStates that Produce Floating Point Error
    df = collect_results(projectdir("Init_States"))

    # Data Frame contains two different rows, each of which produce the same error
    λ = df.λ[2]

    # Boundary Conditions
    BCi = [-0.019488511458668, -0.016033479812051, 0.0,
            8.918881923678198, -4.081793688818725, 0.0,
            1.000000000000000] 

    BCf = [ 0.823385182067467, 0.0, -0.022277556273235,
            0.0, 0.134184170262437, 0.0, 0.0]

    # Time of Flight
    TOF = 8.6404 # [days]

    # Initialize Parameters
    ps = CR3BPIMF_Params("Low Thrust 10", "Earth", "Moon", TOF, BCi, BCf)
    ps.ϵ = 0.0

    # Solve
    sol = nlsolve(only_fj!((F, J, x) -> CR3BPIMF_ShootingFunction_STM!(F, J, x, ps)), λ;
                  method = :newton,
                  xtol = 0.0,
                  ftol = 1e-10,
                  show_trace = true,
                  linesearch = BackTracking())
end

main()