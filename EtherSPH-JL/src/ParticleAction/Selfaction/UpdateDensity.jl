#=
  @ author: bcynuaa
  @ date: 2023-12-05 16:46:24
  @ description:
 =#

function updateDensity!(
    p::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.rho_ += p.drho_ * dt;
    p.drho_ = 0.;
    return nothing;
end