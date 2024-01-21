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
        particles[neighbour.j_].sum_kernel_weight_[Threads.threadid()] += 
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
        particles[neighbour.j_].sum_weighted_scalar_[Threads.threadid()] += 
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_ * 
            getfield(particles[neighbour.i_], value_symbol);
        particles[neighbour.i_].sum_kernel_weight_[Threads.threadid()] += 
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        particles[neighbour.i_].sum_weighted_scalar_[Threads.threadid()] +=
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_ * 
            getfield(particles[neighbour.j_], value_symbol);
    end
    Threads.@threads for particle in particles
        particle.sum_kernel_weight_[1] = sum(particle.sum_kernel_weight_);
        particle.sum_weighted_scalar_[1] = sum(particle.sum_weighted_scalar_);
        particle.sum_kernel_weight_[1] += smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_;
        particle.sum_weighted_scalar_[1] += smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_ * 
            getfield(particle, value_symbol);
        setfield!(particle, value_symbol, particle.sum_weighted_scalar_[1] / particle.sum_kernel_weight_[1]);
        fill!(particle.sum_kernel_weight_, 0.);
        fill!(particle.sum_weighted_scalar_, 0.);
    end
    return nothing;
end

function reconstructVector!(
    particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
    value_symbol::Symbol,
    neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
)::Nothing
    Threads.@threads for neighbour in neighbours
        particles[neighbour.j_].sum_kernel_weight_[Threads.threadid()] += 
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
        particles[neighbour.j_].sum_weighted_vector_[:, Threads.threadid()] .+= 
            neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_ * 
            getfield(particles[neighbour.i_], value_symbol);
        particles[neighbour.i_].sum_kernel_weight_[Threads.threadid()] += 
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        particles[neighbour.i_].sum_weighted_vector_[:, Threads.threadid()] .+=
            neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_ * 
            getfield(particles[neighbour.j_], value_symbol);
    end
    Threads.@threads for particle in particles
        particle.sum_kernel_weight_[1] = sum(particle.sum_kernel_weight_);
        particle.sum_weighted_vector_[:, 1] = sum(particle.sum_weighted_vector_, dims=2);
        particle.sum_kernel_weight_[1] += smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_;
        particle.sum_weighted_vector_[:, 1] .+= smooth_kernel.kernel_0_ * particle.mass_ / particle.rho_ * 
            getfield(particle, value_symbol);
        setfield!(particle, value_symbol, particle.sum_weighted_vector_[:, 1] ./ particle.sum_kernel_weight_[1]);
        fill!(particle.sum_kernel_weight_, 0.);
        fill!(particle.sum_weighted_vector_, 0.);
    end
    return nothing;
end