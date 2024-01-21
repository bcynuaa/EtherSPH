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
    MatrixType <: AbstractMatrix{RealType}
} <: LiquidParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::MatrixType
    rho_::RealType
    drho_::ArrayType
    p_::RealType
    c_::RealType
    mass_::RealType
    gap_::RealType
    sum_kernel_weight_::ArrayType
    sum_weighted_scalar_::ArrayType
    sum_weighted_vector_::MatrixType
end

function LiquidParticleMTh(RealType::DataType, dim::IntType)::LiquidParticleMTh where IntType <: Integer
    n_threads::IntType = Threads.nthreads();
    x_vec = zeros(RealType, dim);
    v_vec = zeros(RealType, dim);
    dv_vec = zeros(RealType, dim, n_threads);
    drho = zeros(RealType, n_threads);
    sum_kernel_weight = zeros(RealType, n_threads);
    sum_weighted_scalar = zeros(RealType, n_threads);
    sum_weighted_vector = zeros(RealType, dim, n_threads);
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