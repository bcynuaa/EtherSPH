#=
  @ author: bcynuaa
  @ date: 2023-12-04 18:55:38
  @ description:
 =#

const kGaussianRadiusRatio = 3.0;
const kGaussianCoefficients = 
[
    1. / sqrt(pi),
    1. / pi,
    1. / sqrt(pi^3)
];

struct GaussianKernel{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::GaussianKernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 3.
        return kernel.sigma_ * exp(-q^2);
    else
        return 0.;
    end
end

function kernelGradient(
    r::RealType,
    kernel::GaussianKernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 3.
        return -2. * kernel.sigma_ / kernel.h_ * q * exp(-q^2);
    else
        return 0.;
    end
end