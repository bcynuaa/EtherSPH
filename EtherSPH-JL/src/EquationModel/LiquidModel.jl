#=
  @ author: bcynuaa
  @ date: 2023-11-24 19:35:16
  @ description:
 =#

abstract type LiquidModel <: AbstractEquationModel end

struct WeakCompressibleLiquidModel <: LiquidModel
    c_0_::Double
    c_02_::Double
    rho_0_::Double
    rho_02_::Double
    p_0_::Double
    mu_::Double
    nu_::Double
    epsilon_rho_::Double
    dim_coefficient_::Double # 2*(dim+2)
    avoid_sigularity_::Double
end

function WeakCompressibleLiquidModel(
  c_0::Double,
  rho_0::Double,
  p_0::Double,
  mu_0::Double,
  dim::Isize
)::WeakCompressibleLiquidModel
    c_02 = c_0 * c_0;
    rho_02 = rho_0 * rho_0;
    nu = mu_0 / rho_0;
    epsilon_rho = 1e-6;
    dim_coefficient = 2. * (dim + 2);
    avoid_sigularity = 1e-2;
    return WeakCompressibleLiquidModel(c_0, c_02, rho_0, rho_02, p_0, mu_0, nu, epsilon_rho, dim_coefficient, avoid_sigularity);
end