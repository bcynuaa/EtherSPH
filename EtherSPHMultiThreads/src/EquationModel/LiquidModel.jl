#=
  @ author: bcynuaa
  @ date: 2024-01-21 15:23:20
  @ description:
 =#

# Weakly Compressible Liquid Model

abstract type WeaklyCompressibleLiquidModel <: LiquidModel end;

struct CommonWeaklyCompressibleLiquidModel{
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
} <: WeaklyCompressibleLiquidModel
    rho_0_::RealType
    c_0_::RealType
    p_0_::RealType
    gamma_::IntType
    b_::RealType # c₀²ρ₀/γ
    mu_0_::RealType
    nu_0_::RealType
    body_force_vec_::ArrayType
end

function CommonWeaklyCompressibleLiquidModel(
    rho_0::RealType, 
    c_0::RealType,
    p_0::RealType,
    gamma::IntType,
    mu_0::RealType,
    body_force_vec::ArrayType
)::CommonWeaklyCompressibleLiquidModel where {
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
}
    b_::RealType = c_0^2 * rho_0 / gamma;
    nu_0_::RealType = mu_0 / rho_0;
    CommonWeaklyCompressibleLiquidModel(
        rho_0, c_0, p_0, 
        gamma, b_, 
        mu_0, nu_0_, 
        body_force_vec
    )
end

function WeaklyCompressibleLiquidModel(
    rho_0::RealType, 
    c_0::RealType,
    p_0::RealType,
    gamma::IntType,
    mu_0::RealType,
    body_force_vec::ArrayType
)::CommonWeaklyCompressibleLiquidModel where {
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
}
    return CommonWeaklyCompressibleLiquidModel(
        rho_0, c_0, p_0, 
        gamma, mu_0, 
        body_force_vec
    )
end

# Weakly Compressible Liquid Thermaldynamic Model

struct ThermalWeaklyCompressibleLiquidModel{
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
} <: WeaklyCompressibleLiquidModel
    rho_0_::RealType
    c_0_::RealType
    p_0_::RealType
    gamma_::IntType
    b_::RealType # c₀²ρ₀/γ
    mu_0_::RealType
    nu_0_::RealType
    body_force_vec_::ArrayType
    t_0_::RealType # reference temperature
    kappa_0_::RealType
    cp_::RealType
    beta_::RealType # ρ-ρ₀ = β(T₀-T)
end

function ThermalWeaklyCompressibleLiquidModel(
    rho_0::RealType, 
    c_0::RealType,
    p_0::RealType,
    gamma::IntType,
    mu_0::RealType,
    body_force_vec::ArrayType,
    t_0::RealType,
    kappa_0::RealType,
    cp::RealType,
    beta::RealType
)::ThermalWeaklyCompressibleLiquidModel where {
    IntType <: Integer, RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType},
}
    b_::RealType = c_0^2 * rho_0 / gamma;
    nu_0_::RealType = mu_0 / rho_0;
    ThermalWeaklyCompressibleLiquidModel(
        rho_0, c_0, p_0, 
        gamma, b_, 
        mu_0, nu_0_, 
        body_force_vec,
        t_0, kappa_0, 
        cp, beta
    )
end