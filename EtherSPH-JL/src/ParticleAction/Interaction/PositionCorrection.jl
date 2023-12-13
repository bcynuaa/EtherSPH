#=
  @ author: bcynuaa
  @ date: 2023-12-05 21:30:20
  @ description:
 =#

function positionCorrection!(
    p_i::ParticleType,
    p_j::ParticleType,
    dt::RealType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    xsph_wc_lm::XSPHWeaklyCompressibleLiquidModel
)::Nothing where {
    RealType <: AbstractFloat,
    ParticleType <: FluidParticle
}
    mean_rho::RealType = (p_i.rho_ + p_j.rho_) / 2.;
    v_correction_vec::typeof(neighbour.v_vec_) = 
        -xsph_wc_lm.epsilon_xsph_ / mean_rho * neighbour.kernel_value_ * neighbour.v_vec_;
    v_i_correction_vec::typeof(v_correction_vec) = p_j.mass_ * v_correction_vec;
    v_j_correction_vec::typeof(v_correction_vec) = -p_i.mass_ * v_correction_vec;
    lock(p_i.lock_) do
        p_i.x_vec_ .+= v_i_correction_vec .* dt;
    end
    lock(p_j.lock_) do
        p_j.x_vec_ .+= v_j_correction_vec .* dt;
    end
    return nothing;
end

function positionCorrection!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    dt::RealType where RealType <: AbstractFloat,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    xsph_wc_lm::XSPHWeaklyCompressibleLiquidModel
)::Nothing
    v_correction_vec::typeof(p_i.v_vec_) = 
        -xsph_wc_lm.epsilon_xsph_ / p_i.rho_ * p_i.mass_ * neighbour.kernel_value_ * neighbour.v_vec_;
    lock(p_i.lock_) do
        p_i.x_vec_ .+= v_correction_vec .* dt;
    end
    return nothing;
end