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
        4 * (mu_i + mu_j) * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
    for (i_dim, dv_i) in enumerate(p_j.mass_ * viscosity_force * neighbour.v_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    for (i_dim, dv_j) in enumerate(-p_i.mass_ * viscosity_force * neighbour.v_vec_)
        Threads.atomic_add!(p_j.dv_vec_[i_dim], dv_j);
    end
    return nothing;
end

function viscosityForce!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    p_j::ParticleType2 where ParticleType2 <: FixedParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing
    sum_rho::typeof(p_i.rho_) = 2 * p_i.rho_;
    mu_i::typeof(sum_rho) = wc_lm.nu_0_ * p_i.rho_;
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        8 * mu_i * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
    for (i_dim, dv_i) in enumerate(p_i.mass_ * viscosity_force * (neighbour.v_vec_ .- p_j.normal_vec_ * dot(neighbour.v_vec_, p_j.normal_vec_)))
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    return nothing;
end

function viscosityForce!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    p_j::ParticleType2 where ParticleType2 <: SolidParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing
    sum_rho::typeof(p_i.rho_) = 2 * p_i.rho_;
    mu_i::typeof(sum_rho) = wc_lm.nu_0_ * p_i.rho_;
    viscosity_force::typeof(sum_rho) = # r ⋅ ∇W = |r|W'
        8 * mu_i * neighbour.r_ * neighbour.kernel_gradient_ / sum_rho^2 / 
            (neighbour.r_^2 + 0.01 * smooth_kernel.h_^2);
    for (i_dim, dv_i) in enumerate(p_i.mass_ * viscosity_force * (neighbour.v_vec_ .- p_j.normal_vec_ * dot(neighbour.v_vec_, p_j.normal_vec_)))
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    for (i_dim, dv_i) in enumerate(p_j.mass_ * viscosity_force * neighbour.v_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    for (i_dim, dv_j) in enumerate(-p_i.mass_ * viscosity_force * neighbour.v_vec_)
        Threads.atomic_add!(p_j.dv_vec_[i_dim], dv_j);
    end
end