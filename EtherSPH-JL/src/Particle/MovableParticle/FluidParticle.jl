#=
  @ author: bcynuaa
  @ date: 2023-11-23 21:54:51
  @ description:
 =#

abstract type FluidParticle <: MovableParticle end

mutable struct LiquidParticle <: FluidParticle
    x_vec_::AbstractVector
    v_vec_::AbstractVector
    dv_vec_::AbstractVector
    rho_::Double
    drho_::Double
    p_::Double
    mass_::Double
end

function LiquidParticle(dim::Isize)::LiquidParticle
    x_vec = zeros(Double, dim);
    v_vec = zeros(Double, dim);
    dv_vec = zeros(Double, dim);
    rho = 0.;
    drho = 0.;
    p = 0.;
    mass = 0.;
    return LiquidParticle(x_vec, v_vec, dv_vec, rho, drho, p, mass);
end