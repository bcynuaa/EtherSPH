#=
  @ author: bcynuaa
  @ date: 2024-01-21 12:18:21
  @ description:
 =#

abstract type FixedParticle <: AbstractParticle end;

abstract type WallParticle <: FixedParticle end; include("./WallParticle.jl");
abstract type VelocityParticle <: FixedParticle end; include("./VelocityParticle.jl");