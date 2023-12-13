#=
  @ author: bcynuaa
  @ date: 2023-12-06 16:11:23
  @ description:
 =#

include("../../src/EtherSPH.jl");
using BenchmarkTools;

p_1 = LiquidParticle(2);
p_1.x_vec_ .= [1., 2.];
p_2 = LiquidParticle(2);
p_2.x_vec_ .= [1.001, 2.002];
p_3 = WallParticle(2);
p_3.x_vec_ .= [1.001, 2.002];

ker = SmoothKernel(0.015, 2, WendlandC2Kernel);

println("="^100);

println("p_1: $p_1");
println("p_2: $p_2");
println("p_3: $p_3");
println("ker: $ker");

println("="^100);

nei1 = CommonNeighbour((1, 2, 0.001414), p_1, p_2, ker);
println("nei CommonNeighbour((1, 2, 0.001414), p_1, p_2, ker): $nei1");
nei2 = CommonNeighbour((1, 3, 0.001414), p_1, p_3, ker);
println("nei CommonNeighbour((1, 3, 0.001414), p_1, p_3, ker): $nei2");

println("="^100);

n = 1500;
pars = [LiquidParticle(2) for i in 1:n];
for i in 1:n
    pars[i].x_vec_ = [rand(), rand()];
end

println("n=$n");
println("findNeighbours(pars, ker)):");
@btime findNeighbours(pars, ker);
println("findNeighbours(pars, pars, ker)):");
@btime findNeighbours(pars, pars, ker);