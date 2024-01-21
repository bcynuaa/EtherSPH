#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:16:12
  @ description:
 =#

 function updateVelocity!(
    particle::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    particle.dv_vec_[:, 1] = sum(particle.dv_vec_, dims=2);
    particle.v_vec_ .+= particle.dv_vec_[:, 1] * dt;
    return nothing;
end

function updateVelocity!(
    particle::ParticleType where ParticleType <: MovableParticle,
    dt::RealType,
    body_force_vec::ArrayType where ArrayType <: AbstractVector{<:RealType}
)::Nothing where RealType <: AbstractFloat
    particle.dv_vec_[:, 1] .+= body_force_vec;
    updateVelocity!(particle, dt);
    return nothing;
end