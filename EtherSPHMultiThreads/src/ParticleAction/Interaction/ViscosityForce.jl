#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:29:35
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
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        4. * (mu_i + mu_j) * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
    p_i.dv_vec_[:, Threads.threadid()] .+= p_j.mass_ * viscosity_force * neighbour.v_vec_;
    p_j.dv_vec_[:, Threads.threadid()] .+= -p_i.mass_ * viscosity_force * neighbour.v_vec_;
    return nothing;
end

function viscosityForce!(
    p_i::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    sum_rho::typeof(p_i.rho_) = 2 * p_i.rho_;
    mu_i::typeof(sum_rho) = wc_lm.nu_0_ * p_i.rho_;
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        4. * 2 * mu_i * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
    p_i.dv_vec_[:, Threads.threadid()] .+= p_i.mass_ * viscosity_force * neighbour.v_vec_;
    return nothing;
end