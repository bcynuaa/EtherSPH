#=
  @ author: bcynuaa
  @ date: 2023-11-27 19:37:40
  @ description:
 =#

function viscosityForce!(
    p_i::AbstractParticle,
    p_j::AbstractParticle,
    particle_neighbour::ParticleNeighbour,
    kernel::AbstractKernel,
    wc_liquid_model::AbstractEquationModel
)::Nothing
    return nothing;
end

function viscosityForce!(
    p_i::FluidParticle,
    p_j::FluidParticle,
    particle_neighbour::ParticleNeighbour,
    kernel::AbstractKernel,
    wc_liquid_model::WeakCompressibleLiquidModel
)::Nothing
    viscosity_force_vec::AbstractVector = 
        wc_liquid_model.dim_coefficient_ * wc_liquid_model.mu_ / wc_liquid_model.rho_02_ .*
            particle_neighbour.kernel_gradient_vec_;
    viscosity_force_vec .*= dot(p_i.v_vec_ .- p_j.v_vec_, particle_neighbour.r_vec_) / 
        (particle_neighbour.r_^2 + wc_liquid_model.avoid_sigularity_ * kernel.h_^2);
    p_i.dv_vec_ .+= p_j.mass_ * viscosity_force_vec;
    p_j.dv_vec_ .-= p_i.mass_ * viscosity_force_vec;
    return nothing;
end

# function viscousityForce!(
#     p_i::FluidParticle,
#     p_j::WallParticle,
#     particle_neighbour::ParticleNeighbour,
#     kernel::AbstractKernel,
#     wc_liquid_model::WeakCompressibleLiquidModel
# )::Nothing
#     viscosity_force_vec::AbstractVector = 
#         wc_liquid_model.dim_coefficient_ * wc_liquid_model.mu_ / wc_liquid_model.rho_02_ .*
#             particle_neighbour.kernel_gradient_vec_;
#     viscosity_force_vec .*= dot(p_i.v_vec_ .- p_j.v_vec_, particle_neighbour.r_vec_) / 
#         (particle_neighbour.r_^2 + wc_liquid_model.avoid_sigularity_ * kernel.h_^2);
#     p_i.dv_vec_ .+= p_i.mass_ * viscosity_force_vec;
#     return nothing;
# end

# function viscousityForce!(
#     p_i::WallParticle,
#     p_j::FluidParticle,
#     particle_neighbour::ParticleNeighbour,
#     kernel::AbstractKernel,
#     wc_liquid_model::WeakCompressibleLiquidModel
# )::Nothing
#     viscosityForce!(p_j, p_i, -particle_neighbour, kernel, wc_liquid_model);
#     return nothing;
# end