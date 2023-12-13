#=
  @ author: bcynuaa
  @ date: 2023-12-06 16:37:37
  @ description:
 =#

include("../../src/EtherSPH.jl");

drfe = DensityReinitializedForwardEuler(0.01, 1000, 100, 30);
@code_warntype DensityReinitializedForwardEuler(0.01, 1000, 100, 30);

isOutputStep(100, drfe)
isDensityReinitializedStep(90, drfe)