#=
  @ author: bcynuaa
  @ date: 2024-02-28 14:19:15
  @ description:
 =#

function thermalConvection!(
    p_i::ParticleType,
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    th_wc_lm::ThermalWeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: FluidParticle
    dq::eltype(p_i.x_vec_) = 1 / th_wc_lm.cp_ * (p_i.kappa_ + p_j.kappa_) * 
        neighbour.r_ * neighbour.kernel_gradient_ / p_i.rho_ / p_j.rho_ /
            (neighbour.r_^2 + 0.01*smooth_kernel.h_^2) * (p_i.t_ - p_j.t_);
    Threads.atomic_add!(p_i.dt_, p_j.mass_ * dq);
    Threads.atomic_add!(p_j.dt_, -p_i.mass_ * dq);
    return nothing;
end

function thermalConvection!(
    p_i::ParticleType1 where ParticleType1 <: FluidParticle,
    p_j::ParticleType2 where ParticleType2 <: ThermostaticWallParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    th_wc_lm::ThermalWeaklyCompressibleLiquidModel
)::Nothing
    dq::eltype(p_i.x_vec_) = 1 / th_wc_lm.cp_ * p_i.kappa_ * 
        neighbour.r_ * neighbour.kernel_gradient_ / p_i.rho_ /
            (neighbour.r_^2 + 0.01*smooth_kernel.h_^2) * (p_i.t_ - p_j.t_);
    Threads.atomic_add!(p_i.dt_, p_i.mass_ * dq);
    return nothing;
end