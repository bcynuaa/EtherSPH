#=
  @ author: bcynuaa
  @ date: 2023-11-23 01:22:25
  @ description:
 =#

const kGaussianRadiusRatio::Double = 3.0;
const kGaussianCoefficients::AbstractVector = 
[
    1. / sqrt(pi),
    1. / pi,
    1. / sqrt(pi^3)
];

struct GaussianKernel <: AbstractKernel
    h_::Double;
    dim_::Isize;
    radius_ratio_::Double;
    influence_radius_::Double;
    sigma_::Double;
end

function GaussianKernel(influence_radius::Double, dim::Isize)::GaussianKernel
    radius_ratio::Double = kGaussianRadiusRatio;
    h::Double = influence_radius / radius_ratio;
    sigma::Double = kGaussianCoefficients[dim] / h^dim;
    return GaussianKernel(h, dim, radius_ratio, influence_radius, sigma);
end

function kernelValue(
    r::Double, 
    kernel::GaussianKernel)::Double
    q::Double = r / kernel.h_;
    if q < 3.0
        return kernel.sigma_ * exp(-q^2);
    else
        return 0.0;
    end
end

function kernelGradient(
    r::Double, 
    kernel::GaussianKernel)::Double
    q::Double = r / kernel.h_;
    if q < 3.0
        return -2.0 * q * kernel.sigma_ * exp(-q^2) / kernel.h_;
    else
        return 0.0;
    end
end