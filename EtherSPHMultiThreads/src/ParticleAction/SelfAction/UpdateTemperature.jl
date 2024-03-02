#=
  @ author: bcynuaa
  @ date: 2024-02-28 15:00:10
  @ description:
 =#

function updateTemperature!(
    p::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat,
)::Nothing
    p.t_ += p.dt_[] * dt;
    p.dt_[] = 0.;
    return nothing;
end