function CR3BPIMF_ShootingFunction_STM!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray}, 
                                        λ::AbstractArray, params::CR3BPIMF_Params)

    ps = deepcopy(params)

    # Time Span
    tspan = (0.0, ps.TOF*3600*24 / ps.eom.TU)
    
    # Initial State Vector
    x0 = Array([view(ps.BCs,:,1); λ; zeros(196)])
    @inbounds for i in 1:14
        x0[14 + (i - 1)*14 + i] = 1.0
    end

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

    # Callback 
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

    if J === nothing
        # Problem Definition
        ff = ODEFunction{true}(CR3BPIMF_EOM!)
        prob = ODEProblem(ff, x0[1:14], tspan, ps, callback = cb)
        alg = Vern9()
    else
        # Problem Definition
        ff = ODEFunction{true}(CR3BPIMF_EOM_W_STM!)
        prob = ODEProblem(ff, x0, tspan, ps, callback = cb)
        alg = TsitPap8()
    end

    # Solve

    sol = try
        solve(prob,
                alg,
                reltol=1e-14, 
                abstol=1e-14, 
                save_everystep=false, 
                save_start=false,
                maxiters=Inf, 
                verbose=true
                ) 
    catch
        safesave(projectdir("Error_Prod_States","states.bson"), Dict(:x0 => x0))
    end

    if !(F === nothing)
        @inbounds for i in 1:6
            F[i] = sol.u[end][i] - params.BCs[i, 2]
        end
        F[7] = sol.u[end][14] - params.BCs[7, 2]

    end
    if !(J === nothing)
        STM = reshape(@view(sol.u[end][15:210]), 14, 14)
        @inbounds for i in 1:6
            J[i,:] .= @view(STM[i,8:14])
        end
        J[7,:] .= @view(STM[14,8:14])
    end
end