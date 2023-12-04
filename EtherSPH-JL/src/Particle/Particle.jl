#=
  @ author: bcynuaa
  @ date: 2023-11-23 21:46:58
  @ description:
 =#

abstract type AbstractParticle end

include("./MovableParticle/MovableParticle.jl");
include("./FixedParticle/FixedParticle.jl");
include("./ParticlePool.jl");