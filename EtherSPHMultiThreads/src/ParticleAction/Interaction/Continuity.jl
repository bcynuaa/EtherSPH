#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:25:45
  @ description:
 =#

function continuity!(
    p_i::ParticleType, 
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour
)::Nothing where ParticleType <: FluidParticle
    drho::typeof(p_i.rho_) = dot(neighbour.v_vec_, neighbour.kernel_gradient_vec_);
    Threads.atomic_add!(p_i.drho_, p_j.mass_ * drho);
    Threads.atomic_add!(p_j.drho_, p_i.mass_ * drho);
    return nothing;
end