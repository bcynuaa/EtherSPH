#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:16:12
  @ description:
 =#

function updateVelocity!(
    particle::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    for i_dim in eachindex(particle.v_vec_)
        particle.v_vec_[i_dim] += particle.dv_vec_[i_dim][] * dt;
    end
    return nothing;
end

function updateVelocity!(
    particle::ParticleType where ParticleType <: MovableParticle,
    dt::RealType,
    body_force_vec::ArrayType where ArrayType <: AbstractVector{<:RealType}
)::Nothing where RealType <: AbstractFloat
    for (i_dim, dv_i) in enumerate(body_force_vec)
        Threads.atomic_add!(particle.dv_vec_[i_dim], dv_i);
    end
    updateVelocity!(particle, dt);
    return nothing;
end