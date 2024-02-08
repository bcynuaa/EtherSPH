#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:34:19
  @ description:
 =#

function wallForce!(
    p_i::ParticleType where ParticleType <: FluidParticle,
    p_j::WallParticleType where WallParticleType <: FixedParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
)::Nothing
    # Roger & Dalrymple, 2008
    psi::eltype(p_i.x_vec_) = abs(dot(neighbour.r_vec_, p_j.normal_vec_));
    xi::typeof(psi) = sqrt(max(neighbour.r_^2 - psi^2, 0.));
    q::typeof(psi) = psi / p_j.gap_;
    if q > 1.
        return nothing;
    end
    if xi > p_j.gap_
        return nothing;
    else
        p_xi::typeof(psi) = abs(0.5 * (1. + cos(pi * xi / p_j.gap_)));
    end
    verticle_u::typeof(psi) = dot(neighbour.v_vec_, p_j.normal_vec_);
    beta = verticle_u > 0. ? 0. : 1.;
    r_psi::typeof(psi) = (0.01 * p_i.c_^2 + beta * p_i.c_ * abs(verticle_u)) / smooth_kernel.h_ * abs(1. - q) / sqrt(q);
    for (i_dim, dv_i) in enumerate(r_psi * p_xi * p_j.normal_vec_)
        Threads.atomic_add!(p_i.dv_vec_[i_dim], dv_i);
    end
    return nothing;
end