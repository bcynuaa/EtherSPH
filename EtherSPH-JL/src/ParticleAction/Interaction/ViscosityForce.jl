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
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        4. * (mu_i + mu_j) * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
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
    dv_vec_i::typeof(p_i.x_vec_) = p_i.mass_ * viscosity_force * neighbour.v_vec_;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    return nothing;
end

function viscosityForce!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    p_j::ParticleType2 where ParticleType2 <: FixedParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    sum_rho::typeof(p_i.rho_) = 2 * p_i.rho_;
    mu_i::typeof(sum_rho) = wc_lm.nu_0_ * p_i.rho_;
    r_vec::typeof(p_i.x_vec_) = p_i.x_vec_ - p_j.x_vec_ - p_j.gap_/2 * p_j.normal_vec_;
    r::eltype(r_vec) = norm(r_vec);
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        4. * 2 * mu_i * r * kernelGradient(r, smooth_kernel) / sum_rho^2 / 
            (r^2 + 0.01 * smooth_kernel.h_^2);
    dv_vec_i::typeof(p_i.x_vec_) = p_i.mass_ * viscosity_force * neighbour.v_vec_;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    return nothing;
end