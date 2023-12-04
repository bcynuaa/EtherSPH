#=
  @ author: bcynuaa
  @ date: 2023-11-28 15:07:11
  @ description:
 =#

abstract type AbstractTimeDiscretization end

function isStepForOutput(
    time_discrete::AbstractTimeDiscretization,
    current_step::Isize
)::Bool
    return current_step % time_discrete.output_step_ == 0;
end

struct ForwardEuler <: AbstractTimeDiscretization
    dt_::Double
    total_step_::Isize
    output_step_::Isize
end