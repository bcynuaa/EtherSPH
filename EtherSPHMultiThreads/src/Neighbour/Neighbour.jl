#=
  @ author: bcynuaa
  @ date: 2024-01-21 16:49:11
  @ description:
 =#

using CellListMap;

abstract type AbstractNeighbour end;

mutable struct CommonNeighbour{
    IntType<:Integer,
    RealType<:AbstractFloat,
    ArrayType<:AbstractVector{RealType}
} <: AbstractNeighbour
    i_::IntType
    j_::IntType
    r_::RealType
    r_vec_::ArrayType
    v_vec_::ArrayType
    kernel_value_::RealType
    kernel_gradient_::RealType
    kernel_gradient_vec_::ArrayType
end

function relativeVelocityVector(
    p_i::ParticleType1 where ParticleType1 <: MovableParticle,
    p_j::ParticleType2 where ParticleType2 <: MovableParticle
)::typeof(p_i.x_vec_)
    return p_i.v_vec_ .- p_j.v_vec_;
end

function relativeVelocityVector(
    p_i::ParticleType1 where ParticleType1 <: MovableParticle,
    p_j::ParticleType2 where ParticleType2 <: FixedParticle
)::typeof(p_i.x_vec_)
    return p_i.v_vec_;
end

function relativeVelocityVector(
    p_i::ParticleType1 where ParticleType1 <: MovableParticle,
    p_j::VelocityParticle
)::typeof(p_i.x_vec_)
    return p_i.v_vec_ .- p_j.v_vec_;
end

function CommonNeighbour(
    neighbour::Tuple{IntType, IntType, RealType},
    p_i::ParticleType1 where ParticleType1 <: AbstractParticle,
    p_j::ParticleType2 where ParticleType2 <: AbstractParticle,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel
)::CommonNeighbour where {
    IntType <: Integer,
    RealType <: AbstractFloat
}   
    i::IntType = neighbour[1];
    j::IntType = neighbour[2];
    r::RealType = neighbour[3];
    r_vec::typeof(p_i.x_vec_) = p_i.x_vec_ .- p_j.x_vec_;
    v_vec::typeof(p_i.x_vec_) = relativeVelocityVector(p_i, p_j);
    kernel_value::typeof(r) = kernelValue(r, smooth_kernel);
    kernel_gradient::typeof(r) = kernelGradient(r, smooth_kernel);
    kernel_gradient_vec::typeof(p_i.x_vec_) = kernel_gradient * normalize(r_vec);
    for i_dim in eachindex(kernel_gradient_vec)
        if isnan(kernel_gradient_vec[i_dim])
            kernel_gradient_vec[i_dim] = 0.;
        end
    end
    return CommonNeighbour(
        i, j, 
        r, r_vec, v_vec,
        kernel_value, kernel_gradient, kernel_gradient_vec
    );
end

function CommonNeighbour(
    neighbour::Tuple{IntType, IntType, RealType},
    p_i::ParticleType1 where ParticleType1 <: AbstractParticle,
    p_j::ParticleType2 where ParticleType2 <: AbstractParticle,
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    common_neighbour::CommonNeighbour
)::Nothing where {
    IntType <: Integer,
    RealType <: AbstractFloat
}   
    common_neighbour.i_ = neighbour[1];
    common_neighbour.j_ = neighbour[2];
    common_neighbour.r_ = neighbour[3];
    common_neighbour.r_vec_ .= p_i.x_vec_ .- p_j.x_vec_;
    common_neighbour.v_vec_ .= relativeVelocityVector(p_i, p_j);
    common_neighbour.kernel_value_ = kernelValue(common_neighbour.r_, smooth_kernel);
    common_neighbour.kernel_gradient_ = kernelGradient(common_neighbour.r_, smooth_kernel);
    common_neighbour.kernel_gradient_vec_ .= common_neighbour.kernel_gradient_ * normalize(common_neighbour.r_vec_);
    for i_dim in eachindex(common_neighbour.kernel_gradient_vec_)
        if isnan(common_neighbour.kernel_gradient_vec_[i_dim])
            common_neighbour.kernel_gradient_vec_[i_dim] = 0.;
        end
    end
    return nothing;
end

function findNeighbours!(
    neighbours::Vector{Tuple{IntType, IntType, RealType}},
    ps_i::ParticleArrayType1 where ParticleArrayType1 <: AbstractVector{<:AbstractParticle},
    ps_j::ParticleArrayType2 where ParticleArrayType2 <: AbstractVector{<:AbstractParticle},
    smooth_kernel::SmoothKernelType where SmoothKernelType <: SmoothKernel,
    neighbours_list::NeighbourArrayType where NeighbourArrayType <: AbstractVector{CommonNeighbour}
)::Nothing where {
    IntType <: Integer,
    RealType <: AbstractFloat
}
    new_neighbours_number::IntType = length(neighbours);
    old_neighbours_number::IntType = length(neighbours_list);
    resize!(neighbours_list, new_neighbours_number);
    @Threads.threads for i_neighbour in eachindex(neighbours)
        if i_neighbour <= old_neighbours_number
            CommonNeighbour(
                neighbours[i_neighbour],
                ps_i[neighbours[i_neighbour][1]],
                ps_j[neighbours[i_neighbour][2]],
                smooth_kernel,
                neighbours_list[i_neighbour]
            );
        else
            neighbours_list[i_neighbour] = CommonNeighbour(
                neighbours[i_neighbour],
                ps_i[neighbours[i_neighbour][1]],
                ps_j[neighbours[i_neighbour][2]],
                smooth_kernel
            );
        end
    end
    return nothing;
end