#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:26:22
  @ description:
 =#

function pressureForce!(
    p_i::ParticleType,
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
)::Nothing where ParticleType <: FluidParticle
    p_rho2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + p_j.p_ / p_j.rho_^2;
    p_rho2 += abs(p_rho2) * 0.01 * neighbour.kernel_value_ / kernelValue((p_i.gap_ + p_j.gap_) / 2., smooth_kernel);
    for (i_dim, dv_i) in enumerate(-p_j.mass_ * p_rho2 * neighbour.kernel_gradient_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    for (i_dim, dv_j) in enumerate(p_i.mass_ * p_rho2 * neighbour.kernel_gradient_vec_)
        Threads.atomic_add!(p_j.dv_vec_[i_dim], dv_j);
    end
    return nothing;
end

function pressureForce!(
    p_i::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    p_rho2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + wc_lm.p_0_ / wc_lm.rho_0_^2;
    p_rho2 += abs(p_rho2) * 0.01 * neighbour.kernel_value_ / kernelValue(p_i.gap_, smooth_kernel);
    for (i_dim, dv_i) in enumerate(-p_i.mass_ * p_rho2 * neighbour.kernel_gradient_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    return nothing;
end

function pressureForce!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    p_j::ParticleType2 where ParticleType2 <: SolidParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModel
)::Nothing
    p_rho2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + wc_lm.p_0_ / wc_lm.rho_0_^2;
    p_rho2 += abs(p_rho2) * 0.01 * neighbour.kernel_value_ / kernelValue((p_i.gap_ + p_j.gap_) / 2., smooth_kernel);
    for (i_dim, dv_i) in enumerate(-p_j.mass_ * p_rho2 * neighbour.kernel_gradient_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    for (i_dim, dv_j) in enumerate(p_i.mass_ * p_rho2 * neighbour.kernel_gradient_vec_)
        Threads.atomic_add!(p_j.dv_vec_[i_dim], dv_j);
    end
    return nothing;
end