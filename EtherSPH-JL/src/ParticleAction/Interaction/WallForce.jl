#=
  @ author: bcynuaa
  @ date: 2023-12-05 18:49:29
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
    if xi < p_j.gap_
        p_xi::typeof(psi) = abs(0.5 * (1. + cos(pi * xi / p_j.gap_))); # taken 2 * gap as frequency
    else
        return nothing;
    end
    verticle_u::typeof(psi) = dot(neighbour.v_vec_, p_j.normal_vec_);
    if verticle_u > 0.
        beta = 0.;
    else
        beta = 1.;
    end
    r_psi::typeof(psi) = (0.01 * p_i.c_^2 + beta * p_i.c_ * abs(verticle_u)) / smooth_kernel.h_ * abs(1. - q) / sqrt(q);
    dv_vec::typeof(p_i.x_vec_) = r_psi * p_xi * p_j.normal_vec_;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec;
    end
    return nothing;
end