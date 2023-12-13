#=
  @ author: bcynuaa
  @ date: 2023-12-05 18:24:02
  @ description:
 =#

function updatePosition!(
    p::ParticleType where ParticleType <: MovableParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.x_vec_ .+= p.v_vec_ .* dt;
    return nothing;
end