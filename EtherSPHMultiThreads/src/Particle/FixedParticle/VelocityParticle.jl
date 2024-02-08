#=
  @ author: bcynuaa
  @ date: 2024-01-21 15:10:17
  @ description:
 =#

mutable struct ConstVelocityParticle{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}
} <: VelocityParticle
    x_vec_::ArrayType
    v_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
end

function ConstVelocityParticle(RealType::DataType, dim::IntType)::ConstVelocityParticle where IntType <: Integer
    ArrayType::DataType = typeof(zeros(RealType, dim));
    x_vec_::ArrayType = zeros(RealType, dim);
    v_vec_::ArrayType = zeros(RealType, dim);
    normal_vec_::ArrayType = zeros(RealType, dim);
    gap_::RealType = 0.;
    return ConstVelocityParticle{RealType, ArrayType}(x_vec_, v_vec_, normal_vec_, gap_);
end