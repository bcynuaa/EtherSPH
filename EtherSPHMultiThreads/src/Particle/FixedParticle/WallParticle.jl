#=
  @ author: bcynuaa
  @ date: 2024-01-21 13:43:25
  @ description:
 =#

mutable struct CompulsiveWallParticle{
    RealType <: AbstractFloat, 
    ArrayType <: AbstractVector{RealType}
} <: WallParticle
    x_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
end

function CompulsiveWallParticle(RealType::DataType, dim::IntType)::CompulsiveWallParticle where IntType <: Integer
    x_vec = zeros(RealType, dim);
    normal_vec = zeros(RealType, dim);
    gap::RealType = 0.;
    CompulsiveWallParticle(x_vec, normal_vec, gap);
end