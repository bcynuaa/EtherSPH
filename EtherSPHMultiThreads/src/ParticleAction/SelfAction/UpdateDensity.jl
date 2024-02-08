#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:07:57
  @ description:
 =#

 function updateDensity!(
    p::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.rho_ += p.drho_[] * dt;
    p.drho_[] = 0.;
    return nothing;
end