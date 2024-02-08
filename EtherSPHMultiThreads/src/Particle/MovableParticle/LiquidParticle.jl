#=
  @ author: bcynuaa
  @ date: 2024-01-21 12:31:59
  @ description:
 =#

"""
# LiquidParticle (multi-thread supported)
MTH -> Multi-Thread
"""
mutable struct LiquidParticleMTh{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}, 
    AtomArrayType <: AbstractVector{Base.Threads.Atomic{RealType}},
} <: LiquidParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::AtomArrayType
    rho_::RealType
    drho_::Base.Threads.Atomic{RealType}
    p_::RealType
    c_::RealType
    mass_::RealType
    gap_::RealType
    sum_kernel_weight_::Base.Threads.Atomic{RealType}
    sum_weighted_scalar_::Base.Threads.Atomic{RealType}
    sum_weighted_vector_::AtomArrayType
end

function LiquidParticleMTh(RealType::DataType, dim::IntType)::LiquidParticleMTh where IntType <: Integer
    x_vec = zeros(RealType, dim);
    v_vec = zeros(RealType, dim);
    dv_vec = [Base.Threads.Atomic{RealType}(0.) for i in 1: dim];
    drho = Base.Threads.Atomic{RealType}(0.);
    sum_kernel_weight = Base.Threads.Atomic{RealType}(0.);
    sum_weighted_scalar = Base.Threads.Atomic{RealType}(0.);
    sum_weighted_vector = [Base.Threads.Atomic{RealType}(0.) for i in 1: dim];
    rho::RealType = 0.;
    p::RealType = 0.;
    c::RealType = 0.;
    mass::RealType = 0.;
    gap::RealType = 0.;
    LiquidParticleMTh(
        x_vec, v_vec, dv_vec, 
        rho, drho, 
        p, c, mass, gap, 
        sum_kernel_weight, 
        sum_weighted_scalar, 
        sum_weighted_vector
    );
end