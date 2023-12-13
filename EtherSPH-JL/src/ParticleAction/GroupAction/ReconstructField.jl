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
    RealType::DataType = typeof(getfield(particles[1], field_symbol));
    weighted_fields = zeros(RealType, length(particles));
    weights = zeros(RealType, length(particles));
    for neighbour in neighbours
        weighted_fields[neighbour.i_] += 
            neighbour.kernel_value_ * getfield(particles[neighbour.j_], field_symbol) * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        weights[neighbour.i_] += neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        weighted_fields[neighbour.j_] += 
            neighbour.kernel_value_ * getfield(particles[neighbour.i_], field_symbol) * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
        weights[neighbour.j_] += neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
    end
    for i_particle in eachindex(particles)
        weighted_fields[i_particle] += 
            smooth_kernel.kernel_0_ * getfield(particles[i_particle], field_symbol) * particles[i_particle].mass_ / particles[i_particle].rho_;
        weights[i_particle] += smooth_kernel.kernel_0_ * particles[i_particle].mass_ / particles[i_particle].rho_;
        setfield!(particles[i_particle], field_symbol, weighted_fields[i_particle] / weights[i_particle]);
    end
    return nothing;
end

function reconstructVector!(
    particles::ParticleArrayType where ParticleArrayType <: AbstractVector{<:AbstractParticle},
    field_symbol::Symbol,
    neighbours::NeighbourArrayType where NeighbourArrayType <: AbstractVector{<:AbstractNeighbour},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
)::Nothing
    RealType::DataType = typeof(getfield(particles[1], field_symbol));
    weighted_fields = zeros(RealType, length(particles), length(getfield(particles[1], field_symbol)));
    weights = zeros(RealType, length(particles));
    for neighbour in neighbours
        weighted_fields[neighbour.i_, :] += 
            neighbour.kernel_value_ * getfield(particles[neighbour.j_], field_symbol) * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        weights[neighbour.i_] += neighbour.kernel_value_ * particles[neighbour.j_].mass_ / particles[neighbour.j_].rho_;
        weighted_fields[neighbour.j_, :] += 
            neighbour.kernel_value_ * getfield(particles[neighbour.i_], field_symbol) * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
        weights[neighbour.j_] += neighbour.kernel_value_ * particles[neighbour.i_].mass_ / particles[neighbour.i_].rho_;
    end
    for i_particle in eachindex(particles)
        weighted_fields[i_particle, :] += 
            smooth_kernel.kernel_0_ * getfield(particles[i_particle], field_symbol) * particles[i_particle].mass_ / particles[i_particle].rho_;
        weights[i_particle] += smooth_kernel.kernel_0_ * particles[i_particle].mass_ / particles[i_particle].rho_;
        setfield!(particles[i_particle], field_symbol, weighted_fields[i_particle, :] / weights[i_particle]);
    end
    return nothing;
end