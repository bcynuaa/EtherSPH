#=
  @ author: bcynuaa
  @ date: 2023-12-05 14:11:09
  @ description:
 =#

abstract type ForwardEuler <: AbstractTimeDiscretization end;

struct DensityReinitializedForwardEuler{IntType <: Integer, RealType <: AbstractFloat} <: ForwardEuler
    dt_::RealType
    total_step_::IntType
    output_step_::IntType
    density_reinitialized_step_::IntType # recommended to be 30 by SPHysics
end

function isDensityReinitializedStep(
    step::IntType where IntType <: Integer,
    time_discretization::DensityReinitializedForwardEuler
)::Bool
    return step % time_discretization.density_reinitialized_step_ == 0;
end