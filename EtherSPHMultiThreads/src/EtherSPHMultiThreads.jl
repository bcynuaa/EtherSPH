module EtherSPHMultiThreads

using LinearAlgebra;
using TypeTree;
using ExportAll;

include("./Kernel/Kernel.jl");
include("./Particle/Particle.jl");
include("./EquationModel/EquationModel.jl");
include("./Neighbour/Neighbour.jl");
include("./ParticleAction/ParticleAction.jl");
include("./DataIO/DataIO.jl");
include("./TimeDiscretization/TimeDiscretization.jl");
include("./Solve/Solve.jl");

function showTypeTree(abstract_type::DataType)::Nothing
    for line::String in tt(abstract_type)
        print(line);
    end
end

println("="^100);

for ether_sph_type in [AbstractKernel, AbstractParticle, AbstractEquationModel, AbstractNeighbour, AbstractTimeDiscretization]
    showTypeTree(ether_sph_type);
    println("-"^80);
end

println("Welcome to EtherSPHMultiThreads.jl!");
println("="^100);

@exportAll;
end # module EtherSPHMultiThreads