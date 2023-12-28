#=
  @ author: bcynuaa
  @ date: 2023-12-05 17:22:50
  @ description:
 =#

function updateVelocity!(
    particle::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    particle.v_vec_ .+= particle.dv_vec_ .* dt;
    return nothing;
end

function updateVelocity!(
    particle::ParticleType where ParticleType <: MovableParticle,
    dt::RealType,
    body_force_vec::ArrayType where ArrayType <: AbstractVector{<:RealType}
)::Nothing where RealType <: AbstractFloat
    particle.dv_vec_ .+= body_force_vec;
    updateVelocity!(particle, dt);
    return nothing;
end