#=
  @ author: bcynuaa
  @ date: 2024-01-21 12:12:02
  @ description:
 =#

abstract type MovableParticle <: AbstractParticle end;

abstract type FluidParticle <: MovableParticle end;
abstract type LiquidParticle <: FluidParticle end; include("./LiquidParticle.jl");
abstract type GasParticle <: FluidParticle end;