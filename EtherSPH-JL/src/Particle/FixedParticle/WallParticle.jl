#=
  @ author: bcynuaa
  @ date: 2023-11-24 00:35:04
  @ description:
 =#

mutable struct WallParticle <: FixedParticle
    x_vec_::AbstractVector
end

function WallParticle(dim::Isize)::WallParticle
    return WallParticle(zeros(Double, dim));
end