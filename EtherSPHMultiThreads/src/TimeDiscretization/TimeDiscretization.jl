#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:49:53
  @ description:
 =#

#=
  @ author: bcynuaa
  @ date: 2023-12-05 14:08:41
  @ description:
 =#

abstract type AbstractTimeDiscretization end;

abstract type FixedTimeIntervalDiscretization <: AbstractTimeDiscretization end; include("./FixedTimeIntervalDiscretization.jl");
abstract type VariableTimeIntervalDiscretization <: AbstractTimeDiscretization end;
    
function isOutputStep(step::IntType, time_discretization::AbstractTimeDiscretization)::Bool where {IntType <: Integer}
    return step % time_discretization.output_step_ == 0;
end

function isChekStep(step::IntType, time_discretization::AbstractTimeDiscretization)::Bool where {IntType <: Integer}
    return step % time_discretization.check_step_ == 0;
end