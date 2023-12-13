#=
  @ author: bcynuaa
  @ date: 2023-12-05 15:46:23
  @ description:
 =#

function viscosityForce!(
    p_i::ParticleType, 
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    sum_rho::typeof(p_i.rho_) = p_i.rho_ + p_j.rho_;
    mu_i::typeof(sum_rho) = wc_lm.nu_0_ * p_i.rho_;
    mu_j::typeof(sum_rho) = wc_lm.nu_0_ * p_j.rho_;
    # r ⋅ ∇W = |r|W'
    viscosity_force::typeof(sum_rho) = 
        4. * (mu_i + mu_j) * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + wc_lm.avoid_singularity_ * smooth_kernel.h_^2);
    # viscosity_force::typeof(sum_rho) = 
    #     4. * (mu_i + mu_j) * neighbour.kernel_gradient_ / sum_rho^2 / neighbour.r_;
    dv_vec_i::typeof(p_i.x_vec_) = p_j.mass_ * viscosity_force * neighbour.v_vec_;
    dv_vec_j::typeof(p_i.x_vec_) = -p_i.mass_ * viscosity_force * neighbour.v_vec_;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    lock(p_j.lock_) do
        p_j.dv_vec_ .+= dv_vec_j;
    end
    return nothing;
end