#=
  @ author: bcynuaa
  @ date: 2023-12-05 14:50:43
  @ description:
 =#

function continuity!(
    p_i::ParticleType, 
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    rho_ij::typeof(p_i.rho_) = p_i.rho_ - p_j.rho_;
    drho::typeof(rho_ij) = dot(neighbour.v_vec_, neighbour.kernel_gradient_vec_);
    # drho_correction::typeof(drho) = 2. * wc_lm.nu_0_ * neighbour.kernel_gradient_ * rho_ij / neighbour.r_;
    drho_correction::typeof(drho) = 0.;
    drho_i::typeof(drho) = p_j.mass_ * (drho + drho_correction);
    drho_j::typeof(drho) = p_i.mass_ * (drho - drho_correction);
    # to avoid data competition, lock here to update the data
    # to speed up, lock's application scope should be as small as possible
    lock(p_i.lock_) do 
        p_i.drho_ += drho_i;
    end
    lock(p_j.lock_) do 
        p_j.drho_ += drho_j;
    end
    return nothing;
end