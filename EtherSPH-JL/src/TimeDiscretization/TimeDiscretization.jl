#=
  @ author: bcynuaa
  @ date: 2023-12-05 14:08:41
  @ description:
 =#

abstract type AbstractTimeDiscretization end;

include("./ForwardEuler.jl");

function isOutputStep(step::IntType, time_discretization::AbstractTimeDiscretization)::Bool where {IntType <: Integer}
    return step % time_discretization.output_step_ == 0;
end