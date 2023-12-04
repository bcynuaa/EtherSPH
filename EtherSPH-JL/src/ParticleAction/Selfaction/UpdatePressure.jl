#=
  @ author: bcynuaa
  @ date: 2023-11-24 16:57:31
  @ description:
 =#

function updatePressure!(p::AbstractParticle, equation_model::AbstractEquationModel)::Nothing
    return nothing;
end

function updatePressure!(p::LiquidParticle, wc_liquid_model::WeakCompressibleLiquidModel)::Nothing
    dp::Double = wc_liquid_model.c_02_ * (p.rho_ - wc_liquid_model.rho_0_);
    p.p_ = wc_liquid_model.p_0_ + dp;
    return nothing;
end