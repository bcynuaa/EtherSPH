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
    p_i.dv_vec_[:, Threads.threadid()] .+= -p_j.mass_ * p_rho2 * neighbour.kernel_gradient_vec_;
    p_j.dv_vec_[:, Threads.threadid()] .+= p_i.mass_ * p_rho2 * neighbour.kernel_gradient_vec_;
    return nothing;
end