#=
  @ author: bcynuaa
  @ date: 2023-12-28 15:39:46
  @ description:
 =#

function saclarKernelCorrection!(
    p_i::ParticleType,
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    scalar_field_symbol::Symbol
)::Nothing where ParticleType <: MovableParticle
    kernel_weight_i::eltype(p_i.x_vec_) = neighbour.kernel_value_ * p_j.mass_ / p_j.rho_;
    kernel_weight_j::typeof(kernel_weight_i) = neighbour.kernel_value_ * p_i.mass_ / p_i.rho_;
    weighted_scalar_i::typeof(kernel_weight_i) = kernel_weight_i * getfield(p_j, scalar_field_symbol);
    weighted_scalar_j::typeof(kernel_weight_i) = kernel_weight_j * getfield(p_i, scalar_field_symbol);
    lock(p_i.lock_) do
        p_i.kernel_weight_ += kernel_weight_i;
        p_i.weighted_scalar_ += weighted_scalar_i;
    end
    lock(p_j.lock_) do
        p_j.kernel_weight_ += kernel_weight_j;
        p_j.weighted_scalar_ += weighted_scalar_j;
    end
    return nothing;
end

function vectorKernelCorrection!(
    p_i::ParticleType,
    p_j::ParticleType,
    neighbour::NeighbourType where NeighbourType <: AbstractNeighbour,
    vector_field_symbol::Symbol
)::Nothing where ParticleType <: MovableParticle
    kernel_weight_i::typeof(p_i.x_vec_) = neighbour.kernel_value_ * p_j.mass_ / p_j.rho_;
    kernel_weight_j::typeof(p_i.x_vec_) = neighbour.kernel_value_ * p_i.mass_ / p_i.rho_;
    weighted_vector_i::typeof(p_i.x_vec_) = kernel_weight_i * getfield(p_j, vector_field_symbol);
    weighted_vector_j::typeof(p_i.x_vec_) = kernel_weight_j * getfield(p_i, vector_field_symbol);
    lock(p_i.lock_) do
        p_i.kernel_weight_ += kernel_weight_i;
        p_i.weighted_vector_ .+= weighted_vector_i;
    end
    lock(p_j.lock_) do
        p_j.kernel_weight_ += kernel_weight_j;
        p_j.weighted_vector_ .+= weighted_vector_j;
    end
    return nothing;
end