#=
  @ author: bcynuaa
  @ date: 2023-12-04 20:39:23
  @ description:
 =#

abstract type FixedParticle <: AbstractParticle end;

include("./WallParticle.jl");
include("./VelocityParticle.jl");