#=
  @ author: bcynuaa
  @ date: 2024-01-03 15:10:05
  @ description:
 =#

mutable struct VelocityParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{RealType}} <: FixedParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
end

function VelocityParticle(dim::IntType)::VelocityParticle where {IntType <: Integer}
    x_vec = zeros(dim);
    v_vec = zeros(dim);
    normal_vec = zeros(dim);
    gap = 0.;
    return VelocityParticle{typeof(gap), typeof(x_vec)}(x_vec, v_vec, normal_vec, gap);
end