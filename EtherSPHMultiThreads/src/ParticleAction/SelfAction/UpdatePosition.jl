#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:18:43
  @ description:
 =#

 function updatePosition!(
    particle::ParticleType where ParticleType <: MovableParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    for i_dim in eachindex(particle.x_vec_)
        particle.x_vec_[i_dim] += particle.v_vec_[i_dim] * dt - particle.dv_vec_[i_dim][] * dt^2 / 2;
        particle.dv_vec_[i_dim][] = 0.;
    end
    return nothing;
end