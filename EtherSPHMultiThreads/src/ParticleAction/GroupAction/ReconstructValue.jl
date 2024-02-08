#=
  @ author: bcynuaa
  @ date: 2024-01-21 17:44:15
  @ description:
 =#

function reconstructScalar!(
    particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
    value_symbol::Symbol,
    neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
)::Nothing
    Threads.@threads for neighbour in neighbours
        Threads.atomic_add!(particles[neighbour.i_].sum_kernel_weight_, 
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_);
        Threads.atomic_add!(particles[neighbour.i_].sum_weighted_scalar_, getfield(particles[neighbour.j_], value_symbol) *
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_);
        Threads.atomic_add!(particles[neighbour.j_].sum_kernel_weight_, 
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_);
        Threads.atomic_add!(particles[neighbour.j_].sum_weighted_scalar_, getfield(particles[neighbour.i_], value_symbol) *
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_);
    end
    Threads.@threads for particle in particles
        Threads.atomic_add!(particle.sum_kernel_weight_, 
            smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_);
        Threads.atomic_add!(particle.sum_weighted_scalar_, getfield(particle, value_symbol) *
            smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_);
        setfield!(particle, value_symbol, particle.sum_weighted_scalar_[] / particle.sum_kernel_weight_[]);
        particle.sum_weighted_scalar_[] = 0.;
        particle.sum_kernel_weight_[] = 0.;
    end
    return nothing;
end