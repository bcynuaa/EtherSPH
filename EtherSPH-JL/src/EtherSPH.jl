#=
  @ author: bcynuaa
  @ date: 2023-11-22 21:39:27
  @ description:
 =#

using LinearAlgebra;
using CellListMap;
using WriteVTK;
using ProgressBars;

const Double = Float64;
const Isize = Int64;

include("./EquationModel/EquationModel.jl");
include("./KernelFunction/KernelFunction.jl");
include("./Particle/Particle.jl");
include("Neighbour/Neighbour.jl");
include("./ParticleAction/ParticleAction.jl");
include("./TimeDiscretization/TimeDiscretization.jl");
include("./DataIO/DataIO.jl");
include("./Solve/Solve.jl");

using TypeTree;

function showTypeTree(some_type::DataType)::Nothing
    for line::String in tt(some_type)
        print(line);
    end
    return nothing;
end

showTypeTree(AbstractKernel);
showTypeTree(AbstractParticle);
showTypeTree(AbstractEquationModel);
showTypeTree(AbstractTimeDiscretization);
showTypeTree(AbstractDataIO);