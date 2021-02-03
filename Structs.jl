# Spacecraft Struct
struct Spacecraft
    # Total Initial Mass
    init_mass::Float64

    # Total Initial Propelant Mass
    prop_mass::Float64

    # Maximum Thrust
    T_max::Float64

    # Specific Impulse
    Isp::Float64

end

# CR3BP Parameters Struct
struct CR3BP_Params
    # Primary and Secondary Body Mass
    m1::Float64
    m2::Float64

    # Radius of Primary and Secondary Bodies
    R1::Float64
    R2::Float64

    # Total Mass
    mtot::Float64

    # Gravitational Parameters
    μ::Float64

    # Distance Between Bodies
    r12::Float64

    # Length Unit
    LU::Float64

    # Time Unit
    TU::Float64

    # Speed Unit
    VU::Float64

    # Mass Unit
    MU::Float64

    # Gravitational Constant
    G::Float64

end

# Constructor
function CR3BP_Params(m1,m2,R1,R2,MU,r12)
    G       = 6.673e-20
    mtot    = m1 + m2
    μ       = m2 / (m1 + m2) 
    LU      = r12
    TU      = 1 / sqrt(G * (m1 + m2) / (r12^3))
    VU      = LU / TU
    MU      = MU
    CR3BP_Params(m1,m2,R1,R2,mtot,μ,r12,LU,TU,VU,MU,G)
end

# CR3BP Indirect Parameters
mutable struct CR3BPIMF_Params
    # Spacecraft Parameters
    sp::Spacecraft

    # EOM Parameters
    eom::CR3BP_Params

    # Homotopic Variable
    ϵ::Float64

    # Thrust Type 
    utype::Int64

    # Time Span
    TOF::Float64

    # Boundry Conditions 
    BCs::AbstractArray{Float64}

    # Max Integration Time
    ITimeMax::Float64

    # Integration T0
    ITime0::Float64
end

# Constructor
function CR3BPIMF_Params(SC::String, PB::String, SB::String,
                         TOF::Number, ICs::AbstractArray, FCs::AbstractArray;
                         ITimeMax::Number = 60.0)

    # Set Spacecraft Parameters 
    if SC == "Orion"
        m0      = 27000.0   # [kg]
        mp      = 8086.4    # [kg]
        T_max   = 26690     # [N]
        Isp     = 315.1     # [s]
        
    elseif SC == "Low Thrust 10"
        m0      = 1500.0    # [kg]
        mp      = 1000.0    # [kg]
        T_max   = 10.0      # [N]
        Isp     = 3000      # [s]
    elseif SC == "Low Thrust 2"
        m0      = 1500.0    # [kg]
        mp      = 1000.0    # [kg]
        T_max   = 2.0      # [N]
        Isp     = 3000      # [s]
    elseif SC == "Low Thrust 1"
        m0      = 1500.0    # [kg]
        mp      = 1000.0    # [kg]
        T_max   = 1.0      # [N]
        Isp     = 3000      # [s]
    elseif SC == "Low Thrust 0.6"
        m0      = 1500.0    # [kg]
        mp      = 1000.0    # [kg]
        T_max   = 0.6       # [N]
        Isp     = 3000      # [s]
    else
        throw(DomainError(SC,"This string does not correspond to an available spacecraft model."))
    end

    # Set CR3BP Parameters
    if PB == "Sun"
        throw(DomainError(PB,"This primary body has not been implemented yet."))
    elseif PB == "Earth"
        m1      = 5.9742e24 # [kg]
        R1      = 6378.137  # [kg]
        
        if SB == "Moon" || SB == "moon"
            m2      = 7.3483e22 # [kg]
            R2      = 1738.0    # [km]
            r12     = 384400    # [km]

        else
            throw(DomainError(SB,"This secondary body has not been implemented for the Earth yet."))
        end
    else
        throw(DomainError(PB,"This primary body has not been implemented yet."))
    end

    return CR3BPIMF_Params(Spacecraft(m0,mp,T_max,Isp),
                           CR3BP_Params(m1,m2,R1,R2,m0,r12),
                           0.0, 0, TOF, SMatrix{length(ICs),2}([ICs FCs]),
                           ITimeMax, time())
end