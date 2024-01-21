#=
  @ author: bcynuaa
  @ date: 2024-01-20 22:42:12
  @ description:
 =#

abstract type SmoothKernel <: AbstractKernel end;

include("./CubicSpline.jl");
include("./Gaussian.jl");
include("./WendlandC2.jl");
include("./WendlandC4.jl");

const kSmoothKernelDict = Dict(
    CubicSpline => (
        2., 
        [
            2.0 / 3.0, 
            10. / 7. / pi, 
            1. / pi
        ]
    ),
    Gaussian => (
        3., 
        [
            1. / sqrt(pi), 
            1. / pi, 
            1. / sqrt(pi^3)
        ]
    ),
    WendlandC2 => (
        2., 
        [
            0.,
            7. / 4. / pi, 
            21. / 16. / pi,
        ]
    ),
    WendlandC4 => (
        2., 
        [
            5. / 8.,
            9. / 4. / pi,
            495. / 256. / pi,
        ]
    ),
);

function SmoothKernel(
    influence_radius::RealType,
    dim::IntType,
    smooth_kernel_type::UnionAll
)::SmoothKernel where {IntType <: Integer, RealType <: AbstractFloat}
    @assert dim == 1 || dim == 2 || dim == 3;
    radius_ratio::RealType = kSmoothKernelDict[smooth_kernel_type][1];
    h::RealType = influence_radius / radius_ratio;
    sigma::RealType = kSmoothKernelDict[smooth_kernel_type][2][dim] / h^dim;
    kernel_0::RealType = sigma;
    return smooth_kernel_type{IntType, RealType}(h, dim, radius_ratio, influence_radius, sigma, kernel_0);
end