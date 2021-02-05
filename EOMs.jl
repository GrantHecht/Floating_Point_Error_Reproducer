# Circular Restricted Three Body Problem EOM with Indirect Min Fuel Control
function CR3BPIMF_EOM!(dx, x, params::CR3BPIMF_Params, t)

    print(t)

    # Compute Requirements
    @fastmath r1      = sqrt((x[1] + params.eom.μ)^2 + x[2]^2 + x[3]^2)
    @fastmath r2      = sqrt((x[1] + params.eom.μ - 1)^2 + x[2]^2 + x[3]^2)
    @fastmath invr13  = r1^(-3)
    @fastmath invr15  = r1^(-5)
    @fastmath invr23  = r2^(-3)
    @fastmath invr25  = r2^(-5)

    # Compute Scalling
    @fastmath T_maxsc = params.sp.T_max * params.eom.TU^2 / (params.eom.MU * params.eom.LU * 1000.0)
    @fastmath c_nsc   = params.sp.Isp * 9.81 * params.eom.TU / (params.eom.LU * 1000.0)

    # Compute Thrust Direction
    #λ_v               = norm(view(x,11:13),2)
    @fastmath λ_v     = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    @fastmath invλ_v  = inv(λ_v) 
    α                 = @SVector [-x[11]*invλ_v, 
                                  -x[12]*invλ_v,
                                  -x[13]*invλ_v]

    # Compute Throttling Factor
    S       = compute_S(x, λ_v, c_nsc)
    u       = compute_u(S, params.utype, params.ϵ)

    # Compute g and h
    g       = @SVector [x[1] - (1 - params.eom.μ)*(x[1] + params.eom.μ)*invr13 - params.eom.μ*(x[1] + params.eom.μ - 1)*invr23,
                        x[2] - (1 - params.eom.μ)*x[2]*invr13 - params.eom.μ*x[2]*invr23,
                        -(1 - params.eom.μ)*x[3]*invr13 - params.eom.μ*x[3]*invr23]
    h       = @SVector [2*x[5], -2*x[4], 0]

    # Compute G and H
    G11     = 1 - (1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*(x[1] + params.eom.μ - 1)^2*invr25;
    G22     = 1 - (1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*x[2]^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*x[2]^2*invr25;
    G33     = -(1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*x[3]^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*x[3]^2*invr25;
    G12     = 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[2]*invr15 + 
                3*params.eom.μ*(x[1] + params.eom.μ - 1)*x[2]*invr25;
    G13     = 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[3]*invr15 + 
                3*params.eom.μ*(x[1] + params.eom.μ - 1)*x[3]*invr25;
    G23     = 3*(1 - params.eom.μ)*x[2]*x[3]*invr15 + 
                3*params.eom.μ*x[2]*x[3]*invr25;

    G       = @SMatrix [ G11 G12 G13;
                         G12 G22 G23;
                         G13 G23 G33]

    H       = @SMatrix [ 0.0 2.0 0.0;
                        -2.0 0.0 0.0;
                         0.0 0.0 0.0]  

    # State/Co-State Derivatives
    dx[7] = -u * T_maxsc / c_nsc
    mag = u * T_maxsc / x[7]
    @inbounds for i in 1:3
        sum1 = 0.0
        sum2 = 0.0
        dx[i] = x[3 + i]
        dx[3 + i] = g[i] + h[i] + α[i]*mag
        for j in 1:3
            sum1 += -G[j, i] * x[10 + j]
            sum2 += -H[j, i] * x[10 + j]
        end
        dx[7 + i] = sum1
        dx[10 + i] = sum2 - x[7 + i]
    end
    dx[14]              = -λ_v*u*T_maxsc / (x[7]^2)
end

# Circular Restricted Three Body Problem EOM with Indirect Min Fuel Control
# and State Transition Matrix
function CR3BPIMF_EOM_W_STM!(dx, x, params::CR3BPIMF_Params, t)

    # Compute Requirements
    @fastmath r1      = sqrt((x[1] + params.eom.μ)^2 + x[2]^2 + x[3]^2)
    @fastmath r2      = sqrt((x[1] + params.eom.μ - 1)^2 + x[2]^2 + x[3]^2)
    @fastmath invr13  = r1^(-3)
    @fastmath invr15  = r1^(-5)
    @fastmath invr23  = r2^(-3)
    @fastmath invr25  = r2^(-5)

    # Compute Scalling
    @fastmath T_maxsc = params.sp.T_max * params.eom.TU^2 / (params.eom.MU*params.eom.LU*1000)
    @fastmath c_nsc   = params.sp.Isp * 9.81 * params.eom.TU / (params.eom.LU * 1000)

    # Compute Thrust Direction
    #λ_v                 = norm(view(x,11:13),2)
    @fastmath λ_v     = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    @fastmath invλ_v  = inv(λ_v) 
    α                 = @SVector [-x[11]*invλ_v, 
                                  -x[12]*invλ_v,
                                  -x[13]*invλ_v]

    # Compute Throttling Factor
    S       = compute_S(x, λ_v, c_nsc)
    u       = compute_u(S, params.utype, params.ϵ)

    # Compute g and h
    g       = @SVector [x[1] - (1 - params.eom.μ)*(x[1] + params.eom.μ)*invr13 - params.eom.μ*(x[1] + params.eom.μ - 1)*invr23,
                        x[2] - (1 - params.eom.μ)*x[2]*invr13 - params.eom.μ*x[2]*invr23,
                        -(1 - params.eom.μ)*x[3]*invr13 - params.eom.μ*x[3]*invr23]
    h       = @SVector [2*x[5], -2*x[4], 0]

    # Compute G and H
    G11     = 1 - (1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*(x[1] + params.eom.μ - 1)^2*invr25;
    G22     = 1 - (1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*x[2]^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*x[2]^2*invr25;
    G33     = -(1 - params.eom.μ)*invr13 + 3*(1 - params.eom.μ)*x[3]^2*invr15 - 
                params.eom.μ*invr23 + 3*params.eom.μ*x[3]^2*invr25;
    G12     = 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[2]*invr15 + 
                3*params.eom.μ*(x[1] + params.eom.μ - 1)*x[2]*invr25;
    G13     = 3*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[3]*invr15 + 
                3*params.eom.μ*(x[1] + params.eom.μ - 1)*x[3]*invr25;
    G23     = 3*(1 - params.eom.μ)*x[2]*x[3]*invr15 + 
                3*params.eom.μ*x[2]*x[3]*invr25;

    G       = @SMatrix [ G11 G12 G13;
                         G12 G22 G23;
                         G13 G23 G33]

    H       = @SMatrix [ 0.0 2.0 0.0;
                        -2.0 0.0 0.0;
                         0.0 0.0 0.0]
        
    # Get STM and Initialize Jacobian
    STM             = reshape(view(x,15:210), 14, 14)
    F               = SizedMatrix{14,14}(zeros(14,14))
    dSTM            = reshape(view(dx,15:210), 14, 14)
    dSTM          .*= 0.0


    # Compute Partials of G 
    invr17  = r1^(-7)
    invr27  = r2^(-7)

    dG11dr          = @SVector [(1 - params.eom.μ)*(params.eom.μ + x[1])*(9*invr15 - 15*(x[1] + params.eom.μ)^2*invr17) +
                    params.eom.μ*(params.eom.μ + x[1] - 1)*(9*invr25 - 15*(x[1] + params.eom.μ - 1)^2*invr27),
                    (1 - params.eom.μ)*x[2]*(3*invr15 - 15*(x[1] + params.eom.μ)^2*invr17) +
                    params.eom.μ*x[2]*(3*invr25 - 15*(x[1] + params.eom.μ - 1)^2*invr27), 
                    (1 - params.eom.μ)*x[3]*(3*invr15 - 15*(x[1] + params.eom.μ)^2*invr17) + 
                    params.eom.μ*x[3]*(3*invr25 - 15*(x[1] + params.eom.μ - 1)^2*invr27)];
                    
    dG12dr          = @SVector [x[2]*(1 - params.eom.μ)*(3*invr15 - 15*(x[1] + params.eom.μ)^2*invr17) + 
                    x[2]*params.eom.μ*(3*invr25 - 15*(x[1] + params.eom.μ - 1)^2*invr27), 
                    (1 - params.eom.μ)*(x[1] + params.eom.μ)*(3*invr15 - 15*x[2]^2*invr17) + 
                    params.eom.μ*(x[1] + params.eom.μ - 1)*(3*invr25 - 15*x[2]^2*invr27), 
                    -15*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[2]*x[3]*invr17 - 
                        15*params.eom.μ*(x[1] + params.eom.μ - 1)*x[2]*x[3]*invr27];

    dG13dr          = @SVector [x[3]*(1 - params.eom.μ)*(3*invr15 - 15*(x[1] + params.eom.μ)^2*invr17) + 
                    x[3]*params.eom.μ*(3*invr25 - 15*(x[1] + params.eom.μ - 1)^2*invr27), 
                    -15*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[2]*x[3]*invr17 - 
                    15*params.eom.μ*(x[1] + params.eom.μ - 1)*x[2]*x[3]*invr27, 
                    (1 - params.eom.μ)*(x[1] + params.eom.μ)*(3*invr15 - 15*x[3]^2*invr17) + 
                    params.eom.μ*(x[1] + params.eom.μ - 1)*(3*invr25 - 15*x[3]^2*invr27)];

    dG22dr          = @SVector [(1 - params.eom.μ)*(x[1] + params.eom.μ)*(3*invr15 - 15*x[2]^2*invr17) + 
                    params.eom.μ*(x[1] + params.eom.μ - 1)*(3*invr25 - 15*x[2]^2*invr27), 
                    (1 - params.eom.μ)*x[2]*(9*invr15 - 15*x[2]^2*invr17) + 
                    params.eom.μ*x[2]*(9*invr25 - 15*x[2]^2*invr27), 
                    (1 - params.eom.μ)*x[3]*(3*invr15 - 15*x[2]^2*invr17) + 
                    params.eom.μ*x[3]*(3*invr25 - 15*x[2]^2*invr27)];

    dG23dr          = @SVector [-15*(1 - params.eom.μ)*(x[1] + params.eom.μ)*x[2]*x[3]*invr17 - 
                    15*params.eom.μ*(x[1] + params.eom.μ - 1)*x[2]*x[3]*invr27, 
                    (1 - params.eom.μ)*x[3]*(3*invr15 - 15*x[2]^2*invr17) + 
                    params.eom.μ*x[3]*(3*invr25 - 15*x[2]^2*invr27), 
                    (1 - params.eom.μ)*x[2]*(3*invr15 - 15*x[3]^2*invr17) + 
                    params.eom.μ*x[2]*(3*invr25 - 15*x[3]^2*invr27)];

    dG33dr          = @SVector [(1 - params.eom.μ)*(x[1] + params.eom.μ)*(3*invr15 - 15*x[3]^2*invr17) + 
                    params.eom.μ*(x[1] + params.eom.μ - 1)*(3*invr25 - 15*x[3]^2*invr27), 
                    (1 - params.eom.μ)*x[2]*(3*invr15 - 15*x[3]^2*invr17) + 
                    params.eom.μ*x[2]*(3*invr25 - 15*x[3]^2*invr27), 
                    (1 - params.eom.μ)*x[3]*(9*invr15 - 15*x[3]^2*invr17) + 
                    params.eom.μ*x[3]*(9*invr25 - 15*x[3]^2*invr27)];
                
    # Compute -dGlamv/dr and place in Jacobian matrix
    @inbounds for i in 1:3
        F[8,i]  = -dG11dr[i]*x[11] - dG12dr[i]*x[12] - dG13dr[i]*x[13]
        F[9,i]  = -dG12dr[i]*x[11] - dG22dr[i]*x[12] - dG23dr[i]*x[13]
        F[10,i] = -dG13dr[i]*x[11] - dG23dr[i]*x[12] - dG33dr[i]*x[13]
    end

    Iden            = SMatrix{3,3}(1.0I)
    Zero            = @SMatrix zeros(3, 3)
    if params.utype != 1
        F[1:3,4:6]              = Iden
        F[4:6,1:3]              = G
        F[4,5]                  = 2.0
        F[5,4]                  = -2.0
        F[8:10,11:13]          .= -transpose(G)
        F[11:13,8:10]          .= -1.0 .* Iden
        F[11,12]                = 2.0
        F[12,11]                = -2.0
        
        compute_G1!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G2!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G7!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G8!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        
    else
        F[1:3,4:6]              = Iden
        F[4:6,1:3]              = G
        F[4,5]                  = 2.0
        F[5,4]                  = -2.0 
        F[8:10,11:13]           = -transpose(G)
        F[11:13,8:10]          .= -1.0 .* Iden
        F[11,12]                = 2.0
        F[12,11]                = -2.0  
        
        compute_G1!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G2!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G3!(F, x, λ_v, T_maxsc, params.ϵ)
        compute_G4!(F, x, λ_v, T_maxsc, params.ϵ)
        compute_G5!(F, x, λ_v, T_maxsc, params.ϵ)
        compute_G6!(F, params.ϵ, c_nsc, T_maxsc)
        compute_G7!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G8!(F, params.utype, x, λ_v, T_maxsc, c_nsc, params.ϵ, u)
        compute_G9!(F, x, λ_v, T_maxsc, params.ϵ)

    end
    # STM Derivatives
    #mul!(dSTM, F, STM, 1.0, 0.0)
    @inbounds Threads.@threads for row in 1:14
        for col in 1:14
            sum = 0.0
            for k in 1:14
                sum += F[row, k] * STM[k, col]
            end
            dSTM[row, col] = sum
        end
    end

    # State/Co-State Derivatives
    #dx[1:3] = @view(x[4:6])
    #dx[4:6] .= g .+ h .+ α.*(u * T_maxsc / x[7])
    dx[7] = -u * T_maxsc / c_nsc
    #@views mul!(dx[8:10],  transpose(G), x[11:13], -1.0, 0)
    #@views mul!(dx[11:13], transpose(H), x[11:13], -1.0, 0) 
    #@views dx[11:13]   .= dx[11:13] .- x[8:10]
    mag = u * T_maxsc / x[7]
    @inbounds for i in 1:3
        sum1 = 0.0
        sum2 = 0.0
        dx[i] = x[3 + i]
        dx[3 + i] = g[i] + h[i] + α[i]*mag
        for j in 1:3
            sum1 += -G[j, i] * x[10 + j]
            sum2 += -H[j, i] * x[10 + j]
        end
        dx[7 + i] = sum1
        dx[10 + i] = sum2 - x[7 + i]
    end
    dx[14]              = -λ_v*u*T_maxsc / (x[7]^2)
end

function compute_S(x::AbstractVector, λ_v, c_nsc)
    return (-λ_v * c_nsc / x[7] - x[14]  + 1.0)
end

function compute_u(S, utype, ϵ)
    if utype == 0
        return 0.0
    elseif utype == 2
        return 1.0
    else
        return (ϵ - S) / (2.0 * ϵ)
    end
end

function compute_G1!(F::AbstractMatrix, utype::Integer, x::AbstractVector, 
                     λ_v, T_maxsc, c_nsc, ϵ, u)
    if utype != 1
        temp = u*T_maxsc / (λ_v*x[7]^2)
        @inbounds for i in 1:3
            F[3 + i, 7] = temp*x[10 + i]
        end

    else
        temp = u*T_maxsc / (λ_v*x[7]^2) + 
               c_nsc*T_maxsc / (2.0*ϵ*x[7]^3)
        @inbounds for i in 1:3
            F[3 + i,7] = temp*x[10 + i]
        end
    end
end

function compute_G2!(F::AbstractMatrix, utype::Integer, x::AbstractVector, 
                     λ_v, T_maxsc, c_nsc, ϵ, u)
    Iden = SMatrix{3,3}(1.0I)
    λ_vinv = 1.0 / λ_v
    if utype != 1
        temp1 = -u*T_maxsc / x[7]
        @inbounds for i in 1:3
            for j in 1:3
                @views F[3+i, 10+j] = temp1*(Iden[i,j]*λ_vinv - x[10 + i]*x[10 + j]*λ_vinv^3)
            end
        end

    else
        temp1                   = -c_nsc*T_maxsc / (2.0*ϵ*x[7]^2)
        temp2                   = -u*T_maxsc / x[7]
        @inbounds for i in 1:3
            for j in 1:3
                F[3+i, 10+j]    = temp1*x[10 + i]*x[10 + j]*λ_vinv^2 + 
                                  temp2*(Iden[i,j]*λ_vinv - x[10 + i]*x[10 + j]*λ_vinv^3)
            end
        end

    end
end

function compute_G3!(F::AbstractMatrix, x::AbstractVector, 
                     λ_v, T_maxsc, ϵ)
    temp = -T_maxsc / (2.0*ϵ*λ_v*x[7])
    @inbounds for i in 1:3
        F[3 + i,14] = temp*x[10 + i]
    end
end

function compute_G4!(F::AbstractMatrix, x::AbstractVector, 
                     λ_v, T_maxsc, ϵ)
    F[7,7] = λ_v*T_maxsc / (2.0*ϵ*x[7]^2)
end

function compute_G5!(F::AbstractMatrix, x::AbstractVector, 
                     λ_v, T_maxsc, ϵ)
    temp = -T_maxsc / (2.0*ϵ*λ_v*x[7])
    @inbounds for i in 1:3
        F[7, 10 + i] = temp*x[10 + i]
    end
end

function compute_G6!(F::AbstractMatrix, ϵ, c_nsc, T_maxsc)
    F[7,14] = -T_maxsc / (2.0*ϵ*c_nsc)
end

function compute_G7!(F::AbstractMatrix, utype::Integer, x::AbstractVector, 
                    λ_v, T_maxsc, c_nsc, ϵ, u)
    if utype != 1
        F[14,7] = 2.0*λ_v*u*T_maxsc / (x[7]^3)
    else
        F[14,7] = 2.0*λ_v*u*T_maxsc / (x[7]^3) + 
                  λ_v^2*c_nsc*T_maxsc / (2.0*ϵ*x[7]^4)
    end
end

function compute_G8!(F::AbstractMatrix, utype::Integer, x::AbstractVector, 
                     λ_v, T_maxsc, c_nsc, ϵ, u)
    if utype != 1
        temp = -u*T_maxsc / (λ_v*x[7]^2)
        @inbounds for i in 1:3
             F[14, 10 + i] = temp*x[10 + i]
        end

    else
        temp = -u*T_maxsc / (λ_v*x[7]^2) - c_nsc*T_maxsc / (2.0*ϵ*x[7]^3)
        @inbounds for i in 1:3
            F[14,10 + i] = temp*x[10 + i]
        end

    end
end

function compute_G9!(F::AbstractMatrix, x::AbstractVector, 
                     λ_v, T_maxsc, ϵ)
    F[14,14] = -λ_v*T_maxsc / (2.0*ϵ*x[7]^2)
end