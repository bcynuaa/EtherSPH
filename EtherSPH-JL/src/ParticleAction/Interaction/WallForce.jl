#=
  @ author: bcynuaa
  @ date: 2023-12-05 18:49:29
  @ description:
 =#

function wallForceEpsilonVerticleDepth(
    p_i::ParticleType where ParticleType <: FluidParticle,
    psi::RealType,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::RealType where RealType <: AbstractFloat
    # this function is needed to be improved
    depth::RealType = p_i.p_ / (wc_lm.rho_0_ * wc_lm.equal_gravity_);
    if depth < 0.
        return 0.02;
    else
        # local_depth::RealType = abs(depth + psi);
        relative_depth::RealType = abs(depth / wc_lm.reference_depth_);
        return min(relative_depth + 0.02, 1.);
    end
end

function wallForceEpsilonVerticleVelocity(
    verticle_u::RealType,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::RealType where RealType <: AbstractFloat
    if verticle_u > 0.
        return 0.;
    else
        epsilon_u::RealType = abs(20. * verticle_u / wc_lm.c_0_);
        return min(epsilon_u, 1.);
    end
end

function wallForce!(
    p_i::ParticleType where ParticleType <: FluidParticle,
    p_j::WallParticle,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    wc_lm::WeaklyCompressibleLiquidModelType where WeaklyCompressibleLiquidModelType <: WeaklyCompressibleLiquidModel
)::Nothing
    #=
    Roger & Dalrymple, 2008
    =#
    psi::eltype(p_i.x_vec_) = abs(dot(neighbour.r_vec_, p_j.normal_vec_));
    xi::typeof(psi) = sqrt(max(neighbour.r_^2 - psi^2, 0.));
    # ! attention here, smooth_kernel.influence_radius_ is not smooth_kernel.h_
    # q::typeof(psi) = psi / smooth_kernel.influence_radius_;
    q::typeof(psi) = psi / p_j.gap_;
    # q::typeof(psi) = psi / (2. * smooth_kernel.h_);
    if q > 1.
        return nothing;
    end
    r_psi::typeof(psi) = 0.01 * p_i.c_^2 * abs(1. - q) / sqrt(q) / (p_j.gap_ / 2.);
    if xi < p_j.gap_
        # taken 2 * gap as frequency
        p_xi = abs(0.5 * (1. + cos(pi * xi / p_j.gap_)));
    else
        p_xi = 0.;
    end
    # p_xi::typeof(psi) = abs(0.5 * (1. + cos(pi * xi / p_j.gap_)));
    verticle_u::typeof(psi) = dot(neighbour.v_vec_, p_j.normal_vec_);
    epsilon_depth::typeof(psi) = wallForceEpsilonVerticleDepth(p_i, psi, wc_lm);
    epsilon_verticle_u::typeof(psi) = wallForceEpsilonVerticleVelocity(verticle_u, wc_lm);
    dv_vec::typeof(p_i.x_vec_) = r_psi * p_xi * (epsilon_depth + epsilon_verticle_u) * p_j.normal_vec_;
    lock(p_i.lock_) do
        p_i.dv_vec_ .+= dv_vec;
    end
    return nothing;
end