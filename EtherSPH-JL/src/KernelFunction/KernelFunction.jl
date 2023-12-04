#=
  @ author: bcynuaa
  @ date: 2023-11-22 21:47:39
  @ description:
 =#

abstract type AbstractKernel end;

include("./CubicSpline.jl");
include("./Gaussian.jl");
include("./WendlandC2.jl");
include("./WendlandC4.jl");

# below 2 functions seem to be useless
# function kernelValue(
#     r_vec::AbstractVector,
#     kernel::AbstractKernel)::Double
#     return kernelValue(norm(r_vec), kernel);
# end

# function kernelGradientVec(
#     r::Double,
#     r_vec::AbstractVector,
#     kernel::AbstractKernel)::AbstractVector
#     return kernelGradient(r, kernel) .* normalize(r_vec);
# end