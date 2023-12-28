#=
  @ author: bcynuaa
  @ date: 2023-12-04 21:33:29
  @ description:
 =#

abstract type LiquidModel <: AbstractEquationModel end;

abstract type WeaklyCompressibleLiquidModel <: LiquidModel end;

struct CommonWeaklyCompressibleLiquidModel{
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
end

function CommonWeaklyCompressibleLiquidModel(
    rho_0::RealType,
    c_0::RealType,
    gamma::RealType,
    mu_0::RealType,
    body_force_vec::ArrayType
)::CommonWeaklyCompressibleLiquidModel where {
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}
}
    b::RealType = c_0^2 * rho_0 / gamma;
    nu_0::RealType = mu_0 / rho_0;
    return CommonWeaklyCompressibleLiquidModel{RealType, ArrayType}(
        rho_0, c_0, gamma, b, 
        mu_0, nu_0, 
        body_force_vec
    );
end