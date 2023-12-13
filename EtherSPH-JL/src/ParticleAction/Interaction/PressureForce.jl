#=
  @ author: bcynuaa
  @ date: 2023-12-05 15:35:39
  @ description:
 =#

function pressureForce!(
    p_i::ParticleType, 
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour
)::Nothing where ParticleType <: FluidParticle
    pressure_force_vec::typeof(p_i.x_vec_) = 
        -(p_i.p_ / p_i.rho_^2 + p_j.p_ / p_j.rho_^2) * neighbour.kernel_gradient_vec_;
    dv_vec_i::typeof(p_i.x_vec_) = p_j.mass_ * pressure_force_vec;
    dv_vec_j::typeof(p_i.x_vec_) = -p_i.mass_ * pressure_force_vec;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    lock(p_j.lock_) do
        p_j.dv_vec_ .+= dv_vec_j;
    end
    return nothing;
end