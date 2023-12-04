#=
  @ author: bcynuaa
  @ date: 2023-11-25 16:51:40
  @ description:
 =#

struct ParticleNeighbour
    i_::Isize
    j_::Isize
    r_::Double
    r_vec_::AbstractVector
    kernel_value_::Double
    kernel_gradient_::Double
    kernel_gradient_vec_::AbstractVector
end

function ParticleNeighbour(
    i::Isize,
    j::Isize,
    p_i::AbstractParticle,
    p_j::AbstractParticle,
    kernel::AbstractKernel
)::ParticleNeighbour
    r_vec::AbstractVector = p_i.x_vec_ .- p_j.x_vec_;
    r::Double = norm(r_vec);
    kernel_value::Double = kernelValue(r, kernel);
    kernel_gradient::Double = kernelGradient(r, kernel);
    kernel_gradient_vec::AbstractVector = kernel_gradient * normalize(r_vec);
    return ParticleNeighbour(i, j, r, r_vec, kernel_value, kernel_gradient, kernel_gradient_vec);
end

function ParticleNeighbour(
    i::Isize,
    j::Isize,
    p_i::AbstractParticle,
    p_j::AbstractParticle,
    r::Double,
    kernel::AbstractKernel
)::ParticleNeighbour
    r_vec::AbstractVector = p_i.x_vec_ .- p_j.x_vec_;
    kernel_value::Double = kernelValue(r, kernel);
    kernel_gradient::Double = kernelGradient(r, kernel);
    kernel_gradient_vec::AbstractVector = kernel_gradient * normalize(r_vec);
    return ParticleNeighbour(i, j, r, r_vec, kernel_value, kernel_gradient, kernel_gradient_vec);
end

function ParticleNeighbour(
    neighbor::Tuple{Isize, Isize, Double},
    particle_pool::ParticlePool,
    kernel::AbstractKernel
)::ParticleNeighbour
    return ParticleNeighbour(
        neighbor[1],
        neighbor[2],
        particle_pool.particles_[neighbor[1]],
        particle_pool.particles_[neighbor[2]],
        neighbor[3],
        kernel
    );
end

function Base.:-(
    particle_neighbour::ParticleNeighbour,
)::ParticleNeighbour
    return ParticleNeighbour(
        particle_neighbour.j_,
        particle_neighbour.i_,
        particle_neighbour.r_,
        -particle_neighbour.r_vec_,
        particle_neighbour.kernel_value_,
        particle_neighbour.kernel_gradient_,
        -particle_neighbour.kernel_gradient_vec_
    );
end

function findNeighbours(
    particle_pool::ParticlePool,
    kernel::AbstractKernel
)::AbstractVector{ParticleNeighbour}
    neighbours::AbstractVector{Tuple{Isize, Isize, Double}} = 
        neighborlist([particle.x_vec_ for particle in particle_pool.particles_], kernel.influence_radius_);
    return [ParticleNeighbour(neighbor, particle_pool, kernel) for neighbor in neighbours];
end