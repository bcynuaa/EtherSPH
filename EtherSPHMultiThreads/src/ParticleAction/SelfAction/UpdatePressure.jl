#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:11:13
  @ description:
 =#

 function updatePressure!(
    p::ParticleType where ParticleType <: FluidParticle,
    wc_lm::WeaklyCompressibleLiquidModel
)::Nothing
    #=
    p = c₀²ρ₀/γ * [(ρ/ρ₀)^γ - 1]
    c = √(dp/dρ) = c₀ * (ρ/ρ₀)^(γ-1)/2
    =#
    relative_rho::typeof(p.rho_) = p.rho_ / wc_lm.rho_0_;
    p.p_ = wc_lm.b_ * (relative_rho^wc_lm.gamma_ - 1.) + wc_lm.p_0_;
    p.c_ = wc_lm.c_0_ * sqrt(relative_rho^(wc_lm.gamma_ - 1));
    return nothing;
end