#=
  @ author: bcynuaa
  @ date: 2023-12-04 21:33:29
  @ description:
 =#

abstract type LiquidModel <: AbstractEquationModel end;

abstract type WeaklyCompressibleLiquidModel <: LiquidModel end;

struct XSPHWeaklyCompressibleLiquidModel{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}
} <: WeaklyCompressibleLiquidModel
    rho_0_::RealType
    c_0_::RealType
    gamma_::RealType
    b_::RealType # c₀²ρ₀/γ
    mu_0_::RealType
    nu_0_::RealType
    body_force_vec_::ArrayType
    equal_gravity_::RealType
    epsilon_xsph_::RealType # modification while moving particles
    reference_depth_::RealType
    avoid_singularity_::RealType
end

function XSPHWeaklyCompressibleLiquidModel(
    rho_0::RealType,
    c_0::RealType,
    gamma::RealType,
    mu_0::RealType,
    body_force_vec::ArrayType,
    reference_depth::RealType
)::XSPHWeaklyCompressibleLiquidModel where {
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}
}
    b::RealType = c_0^2 * rho_0 / gamma;
    nu_0::RealType = mu_0 / rho_0;
    equal_gravity::RealType = norm(body_force_vec);
    epsilon_xsph::RealType = 0.5;
    avoid_singularity::RealType = 1e-2;
    return XSPHWeaklyCompressibleLiquidModel{RealType, ArrayType}(
        rho_0, c_0, gamma, b, 
        mu_0, nu_0, 
        body_force_vec, equal_gravity,
        epsilon_xsph, reference_depth, 
        avoid_singularity
    );
end