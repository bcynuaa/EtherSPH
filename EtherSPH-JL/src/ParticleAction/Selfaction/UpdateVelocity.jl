#=
  @ author: bcynuaa
  @ date: 2023-12-05 17:22:50
  @ description:
 =#

function updateVelocity!(
    p::ParticleType where ParticleType <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat
)::Nothing
    p.v_vec_ .+= p.dv_vec_ .* dt;
    p.dv_vec_ .= 0.;
    return nothing;
end

function updateVelocity!(
    p::ParticleType where ParticleType <: MovableParticle,
    dt::RealType,
    body_force_vec::ArrayType where ArrayType <: AbstractVector{<:RealType}
)::Nothing where RealType <: AbstractFloat
    p.dv_vec_ .+= body_force_vec;
    updateVelocity!(p, dt);
    return nothing;
end