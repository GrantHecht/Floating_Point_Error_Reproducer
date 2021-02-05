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

function main()

    # Load Error Producing State Vectors 
    # (Each row of data frame contains an error producing state vector)
    df = collect_results(projectdir("Error_Prod_States"))
    x0 = df.x0[2][1:14]

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

    # Time Span
    tspan = (0.0, TOF*3600*24 / ps.eom.TU)

    # Set Utype
    c_nsc = ps.sp.Isp * 9.81 * ps.eom.TU / (ps.eom.LU * 1000)
    λ_v = norm(view(x0,11:13),2)
    S = compute_S(x0, λ_v, c_nsc)
    if S > ps.ϵ
        ps.utype = 0
    elseif S < -ps.ϵ
        ps.utype = 2
    else
        ps.utype = 1
    end

     # Callbacks
     cb1 = ContinuousCallback(CR3BPIMF_Switching_Condition,
                            CR3BPIMF_Switching_Affect!;
                            idxs = 1:14,
                            rootfind = true,
                            affect_neg! = CR3BPIMF_Switching_Affect!,
                            interp_points = 50,
                            abstol = 10eps(),
                            reltol = 0)

    cb2 = AutoAbstol(false)
    cb = CallbackSet(cb1,cb2)

    # Problem Definition
    ff = ODEFunction{true}(CR3BPIMF_EOM!)
    prob = ODEProblem(ff, x0[1:14], tspan, ps, callback = cb)
    alg = Vern9()

    # Solve
    sol = solve(prob,
                alg,
                reltol=1e-14, 
                abstol=1e-14, 
                save_everystep=false, 
                save_start=false,
                maxiters=Inf, 
                verbose=true
                ) 

end

main()