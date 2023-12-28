#=
  @ author: bcynuaa
  @ date: 2023-12-05 14:50:43
  @ description:
 =#

function continuity!(
    p_i::ParticleType, 
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour
)::Nothing where ParticleType <: FluidParticle
    rho_ij::typeof(p_i.rho_) = p_i.rho_ - p_j.rho_;
    drho::typeof(rho_ij) = dot(neighbour.v_vec_, neighbour.kernel_gradient_vec_);
    drho_i::typeof(drho) = p_j.mass_ * drho;
    drho_j::typeof(drho) = p_i.mass_ * drho;
    lock(p_i.lock_) do 
        p_i.drho_ += drho_i;
    end
    lock(p_j.lock_) do 
        p_j.drho_ += drho_j;
    end
    return nothing;
end