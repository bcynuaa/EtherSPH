#=
  @ author: bcynuaa
  @ date: 2023-12-04 20:37:49
  @ description:
 =#

abstract type AbstractParticle end;

include("./MovableParticle/MovableParticle.jl");
include("./FixedParticle/FixedParticle.jl");