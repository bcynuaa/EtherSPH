#=
  @ author: bcynuaa
  @ date: 2023-12-04 20:40:53
  @ description:
 =#

mutable struct WallParticle{RealType <: AbstractFloat, ArrayType <: AbstractVector{RealType}} <: FixedParticle
    x_vec_::ArrayType
    normal_vec_::ArrayType
    gap_::RealType
end

function WallParticle(dim::IntType)::WallParticle where {IntType <: Integer}
    x_vec = zeros(dim);
    normal_vec = zeros(dim);
    gap = 0.;
    return WallParticle{typeof(gap), typeof(x_vec)}(x_vec, normal_vec, gap);
end