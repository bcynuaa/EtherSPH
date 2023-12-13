#=
  @ author: bcynuaa
  @ date: 2023-12-04 17:14:36
  @ description:
 =#

using LinearAlgebra;

include("./Kernel/Kernel.jl");
include("./Particle/Particle.jl");
include("./EquationModel/EquationModel.jl");
include("./Neighbour/Neighbour.jl");
include("./ParticleAction/ParticleAction.jl");
include("./TimeDiscretization/TimeDiscretization.jl");
include("./DataIO/DataIO.jl");
include("./Solve/Solve.jl");

using TypeTree;

function showTypeTree(abstract_type::DataType)::Nothing
    for line::String in tt(abstract_type)
        print(line);
    end
end

showTypeTree(AbstractKernel);
showTypeTree(AbstractParticle);
showTypeTree(AbstractEquationModel);
showTypeTree(AbstractNeighbour);
showTypeTree(AbstractTimeDiscretization);