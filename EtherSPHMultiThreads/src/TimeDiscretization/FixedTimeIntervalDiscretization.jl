#=
  @ author: bcynuaa
  @ date: 2024-01-21 18:52:53
  @ description:
 =#

struct FixedDensityReinitializedForwardEuler{IntType <: Integer, RealType <: AbstractFloat} <: FixedTimeIntervalDiscretization
    dt_::RealType
    total_step_::IntType
    output_step_::IntType
    density_reinitialized_step_::IntType # recommended to be 30 by SPHysics
    check_step_::IntType
end

function isDensityReinitializedStep(
    step::IntType where IntType <: Integer,
    time_discretization::FixedDensityReinitializedForwardEuler
)::Bool
    return step % time_discretization.density_reinitialized_step_ == 0;
end