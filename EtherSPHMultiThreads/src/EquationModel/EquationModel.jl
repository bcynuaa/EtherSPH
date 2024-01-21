#=
  @ author: bcynuaa
  @ date: 2024-01-21 15:21:51
  @ description:
 =#

abstract type AbstractEquationModel end;

abstract type LiquidModel <: AbstractEquationModel end; include("./LiquidModel.jl");