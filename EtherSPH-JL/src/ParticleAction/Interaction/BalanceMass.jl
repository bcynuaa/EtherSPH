#=
  @ author: bcynuaa
  @ date: 2023-11-24 21:34:11
  @ description:
 =#

function balanceMass!(
    p_i::AbstractParticle,
    p_j::AbstractParticle,
    particle_neighbour::ParticleNeighbour,
    equation_model::AbstractEquationModel
)::Nothing
    return nothing;
end

#=
  `` \frac{\partial \rho_i}{\partial t} = - \sum_j m_j \frac{\partial W_{ij}}{\partial r_{ij}} \cdot \vec{v}_{ij} ``
=#
function balanceMass!(
    p_i::FluidParticle,
    p_j::FluidParticle,
    particle_neighbour::ParticleNeighbour,
    wc_liquid_model::WeakCompressibleLiquidModel
)::Nothing
    rho_ij::Double = p_i.rho_ - p_j.rho_;
    v_ij_vec::AbstractVector = p_i.v_vec_ .- p_j.v_vec_;
    drho::Double = dot(v_ij_vec, particle_neighbour.kernel_gradient_vec_);
    drho_correction::Double = 
        2. * wc_liquid_model.epsilon_rho_ * particle_neighbour.kernel_gradient_ * rho_ij / particle_neighbour.r_;
    p_i.drho_ += p_j.mass_ * (drho + drho_correction);
    p_j.drho_ += p_i.mass_ * (drho - drho_correction);
    return nothing;
end