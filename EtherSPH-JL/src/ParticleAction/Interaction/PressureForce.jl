#=
  @ author: bcynuaa
  @ date: 2023-12-05 15:35:39
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
    pressure_force_vec::typeof(p_i.x_vec_) = -p_rho2 * neighbour.kernel_gradient_vec_;
    dv_vec_i::typeof(p_i.x_vec_) = p_j.mass_ * pressure_force_vec;
    dv_vec_j::typeof(p_i.x_vec_) = -p_i.mass_ * pressure_force_vec;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    lock(p_j.lock_) do
        p_j.dv_vec_ .+= dv_vec_j;
    end
    return nothing;
end

function pressureForce!(
    p_i::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing where ParticleType <: MovableParticle
    p_rho2::typeof(p_i.p_) = p_i.p_ / p_i.rho_^2 + wc_lm.p_0_ / wc_lm.rho_0_^2;
    # p_rho2 += abs(p_rho2) * 0.01 * neighbour.kernel_value_ / kernelValue(p_i.gap_, smooth_kernel);
    pressure_force_vec::typeof(p_i.x_vec_) = -p_rho2 * neighbour.kernel_gradient_vec_;
    dv_vec_i::typeof(p_i.x_vec_) = p_i.mass_ * pressure_force_vec;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec_i;
    end
    return nothing;
end