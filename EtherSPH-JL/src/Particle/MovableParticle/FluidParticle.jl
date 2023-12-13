#=
  @ author: bcynuaa
  @ date: 2023-12-04 20:49:07
  @ description:
 =#

abstract type FluidParticle <: MovableParticle end;

mutable struct LiquidParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{RealType}} <: FluidParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    dv_vec_::ArrayType
    rho_::RealType
    drho_::RealType
    p_::RealType
    c_::RealType
    mass_::RealType
    lock_::ReentrantLock
end

function LiquidParticle(dim::IntType)::LiquidParticle where {IntType <: Integer}
    x_vec = zeros(dim);
    v_vec = zeros(dim);
    dv_vec = zeros(dim);
    rho = 0.;
    drho = 0.;
    p = 0.;
    c = 0.;
    mass = 0.;
    lock::ReentrantLock = Threads.ReentrantLock();
    return LiquidParticle{typeof(mass), typeof(x_vec)}(x_vec, v_vec, dv_vec, rho, drho, p, c, mass, lock);
end