function CR3BPIMF_Switching_Condition(x, t, integrator)

    # Compute Requirements
    c_nsc = integrator.p.sp.Isp * 9.81 * integrator.p.eom.TU / (integrator.p.eom.LU * 1000)
    λ_v = norm(view(x,11:13),2)
    S = compute_S(x, λ_v, c_nsc)

    # Switching Condition
    if integrator.p.ϵ != 0.0
        if integrator.p.utype == 0
            return (S - integrator.p.ϵ)
        elseif integrator.p.utype == 2
            return (S + integrator.p.ϵ)
        else
            return (abs(S) - integrator.p.ϵ)
        end
    else
        return S
    end
end

function CR3BPIMF_Switching_Affect!(integrator)

    # Switching Affect
    if integrator.p.ϵ != 0.0
        if integrator.p.utype != 1
            integrator.p.utype = 1
        else
            # Compute Switching Function
            @fastmath c_nsc = integrator.p.sp.Isp * 9.81 * integrator.p.eom.TU / (integrator.p.eom.LU * 1000)
            @fastmath λ_v = norm(view(integrator.u,11:13),2)
            @fastmath S = compute_S(integrator.u, λ_v, c_nsc)

            if S < 0.0
                integrator.p.utype = 2
            else 
                integrator.p.utype = 0
            end
        end
    else
        flag = false

        # If Propagating STM
        if length(integrator.u) == 210

            # Initialize Requirements
            dy⁻ = MVector{14}(zeros(14))
            dy⁺ = MVector{14}(zeros(14))
            flag = true
            
            # Compute dy⁻
            CR3BPIMF_EOM!(dy⁻, MVector{14}(@views(integrator.u[1:14])), integrator.p, 0.0)

            # Compute Requirements
            @fastmath c_nsc = integrator.p.sp.Isp * 9.81 * integrator.p.eom.TU / (integrator.p.eom.LU * 1000)
            @fastmath λ_v = norm(view(integrator.u,11:13),2)

            # Compute ∂S and dSinv
            ∂S = MVector{14}([zeros(6); 
                                λ_v*c_nsc/(integrator.u[7]^2);
                                zeros(3);
                                @views(integrator.u[11:13]).*(-c_nsc/(λ_v*integrator.u[7]));
                                -1.0])
            dSinv = (integrator.u[8] - 2*integrator.u[12])*integrator.u[11] +
                    (integrator.u[9] + 2*integrator.u[11])*integrator.u[12] +
                    (integrator.u[10])*integrator.u[13]
            dSinv = λ_v*integrator.u[7]/(c_nsc*dSinv)
            
        end

        # Set Utype
        if integrator.p.utype == 0
            integrator.p.utype = 2
        else
            integrator.p.utype = 0
        end

        # If Propagating STM
        if flag

            # Compute dy⁺
            CR3BPIMF_EOM!(dy⁺, MVector{14}(@views(integrator.u[1:14])), integrator.p, 0.0)

            # Compute Φ and Ψ
            uSTM = reshape(view(integrator.u,15:210), 14, 14)
            Φ = SizedMatrix{14,14}(zeros(14,14))
            Ψ = SizedMatrix{14,14}(zeros(14,14))
            Threads.@threads for row in 1:14
                for col in 1:14
                    Φ[row, col] = uSTM[row, col]
                    Ψ[row, col] = I[row,col] + (dy⁺[row] - dy⁻[row]) * ∂S[col]*dSinv
                end
            end
            #Ψ .= SizedMatrix{14,14}(I) .+ (dy⁺ .- dy⁻)*transpose(∂S.*dSinv)

            # Multiply Ψ and Φ
            Threads.@threads for row in 1:14
                for col in 1:14
                    sum = 0.0
                    for k in 1:14
                        sum += Ψ[row, k] * Φ[k, col]
                    end
                    uSTM[row, col] = sum
                end
            end

        end
    end

    # Reset dt
    auto_dt_reset!(integrator)
    #set_proposed_dt!(integrator, get_proposed_dt(integrator) / 100.0)
end