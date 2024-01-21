#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:18:43
  @ description:
 =#

 function updatePosition!(
    p::ParticleType where ParticleType <: MovableParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.x_vec_ .+= p.v_vec_ .* dt .- p.dv_vec_[:, 1] .* dt^2 / 2.;
    fill!(p.dv_vec_, 0.);
    return nothing;
end