#=
  @ author: bcynuaa
  @ date: 2023-11-25 16:43:33
  @ description:
 =#

mutable struct ParticlePool
    particles_::AbstractVector{AbstractParticle};
    type_index_dict_::Dict{DataType, AbstractVector};
end

function ParticlePool(
    particles::AbstractVector{AbstractParticle}
)::ParticlePool
    type_index_dict::Dict{DataType, AbstractVector} = Dict();
    n_particles::Isize = length(particles);
    for i = 1:n_particles
        particle_type::DataType = typeof(particles[i]);
        if haskey(type_index_dict, particle_type)
            push!(type_index_dict[particle_type], i);
        else
            type_index_dict[particle_type] = [i];
        end
    end
    return ParticlePool(particles, type_index_dict);
end