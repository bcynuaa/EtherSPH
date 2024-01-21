#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:07:57
  @ description:
 =#

 function updateDensity!(
    p::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.drho_[1] = sum(p.drho_);
    p.rho_ += p.drho_[1] * dt;
    fill!(p.drho_, 0.);
    return nothing;
end