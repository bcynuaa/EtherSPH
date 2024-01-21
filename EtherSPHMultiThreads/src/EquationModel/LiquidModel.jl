#=
  @ author: bcynuaa
  @ date: 2024-01-21 15:23:20
  @ description:
 =#

struct WeaklyCompressibleLiquidModel{
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
} <: LiquidModel
    rho_0_::RealType
    c_0_::RealType
    p_0_::RealType
    gamma_::IntType
    b_::RealType # c₀²ρ₀/γ
    mu_0_::RealType
    nu_0_::RealType
    body_force_vec_::ArrayType
end

function WeaklyCompressibleLiquidModel(
    rho_0::RealType, 
    c_0::RealType,
    p_0::RealType,
    gamma::IntType,
    mu_0::RealType,
    body_force_vec::ArrayType
)::WeaklyCompressibleLiquidModel where {
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
}
    b_::RealType = c_0^2 * rho_0 / gamma;
    nu_0_::RealType = mu_0 / rho_0;
    WeaklyCompressibleLiquidModel(
        rho_0, c_0, p_0, 
        gamma, b_, 
        mu_0, nu_0_, 
        body_force_vec
    )
end