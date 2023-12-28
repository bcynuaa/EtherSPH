#=
  @ author: bcynuaa
  @ date: 2023-12-09 01:38:38
  @ description:
 =#

function reconstructScalar!(
    particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
    field_symbol::Symbol,
    neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
)::Nothing
    Threads.@threads for neighbour in neighbours
        saclarKernelCorrection!(
            particles[neighbour.i_],
            particles[neighbour.j_],
            neighbour,
            field_symbol
        );
    end
    Threads.@threads for i_particle in eachindex(particles)
        particles[i_particle].kernel_weight_ += smooth_kernel.kernel_0_ * particles[i_particle].mass_ / particles[i_particle].rho_;
        particles[i_particle].weighted_scalar_ += 
            smooth_kernel.kernel_0_ * getfield(particles[i_particle], field_symbol) * 
                particles[i_particle].mass_ / particles[i_particle].rho_;
        setfield!(particles[i_particle], field_symbol, particles[i_particle].weighted_scalar_ / particles[i_particle].kernel_weight_);
        particles[i_particle].kernel_weight_ = 0.;
        particles[i_particle].weighted_scalar_ = 0.;
    end
    return nothing;
end

function reconstructVector!(
    particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
    field_symbol::Symbol,
    neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
)::Nothing
    Threads.@threads for neighbour in neighbours
        vectorKernelCorrection!(
            particles[neighbour.i_],
            particles[neighbour.j_],
            neighbour,
            field_symbol
        );
    end
    Threads.@threads for i_particle in eachindex(particles)
        particles[i_particle].kernel_weight_ += smooth_kernel.kernel_0_ * particles[i_particle].mass_ / particles[i_particle].rho_;
        particles[i_particle].weighted_vector_ .+= 
            smooth_kernel.kernel_0_ * getfield(particles[i_particle], field_symbol) * 
                particles[i_particle].mass_ / particles[i_particle].rho_;
        setfield!(particles[i_particle], field_symbol, particles[i_particle].weighted_vector_ / particles[i_particle].kernel_weight_);
        particles[i_particle].kernel_weight_ = 0.;
        particles[i_particle].weighted_vector_ .= 0.;
    end
    return nothing;
end