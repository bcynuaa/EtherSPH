#=
  @ author: bcynuaa
  @ date: 2023-12-06 16:35:53
  @ description:
 =#

include("../../src/EtherSPH.jl");
using BenchmarkTools;

rho_0 = 1000.;
c_0 = 1500.;
gamma = 7.;
mu_0 = 8.9e-4;
body_force_vec = [0., -9.81];
xsph_wc_lm = XSPHWeaklyCompressibleLiquidModel(rho_0, c_0, gamma, mu_0, body_force_vec);
@code_warntype XSPHWeaklyCompressibleLiquidModel(rho_0, c_0, gamma, mu_0, body_force_vec);