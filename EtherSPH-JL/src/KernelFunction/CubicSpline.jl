#=
  @ author: bcynuaa
  @ date: 2023-11-22 21:49:54
  @ description:
 =#

const kCubicSplineRadiusRatio::Double = 2.0;
const kCubicSplineCoefficients::AbstractVector = 
[
    2.0 / 3.0,
    10. / 7. / pi,
    1. / pi,
];

struct CubicSplineKernel <: AbstractKernel
    h_::Double;
    dim_::Isize;
    radius_ratio_::Double;
    influence_radius_::Double;
    sigma_::Double;
end

function CubicSplineKernel(influence_radius::Double, dim::Isize)::CubicSplineKernel
    radius_ratio::Double = kCubicSplineRadiusRatio;
    h::Double = influence_radius / radius_ratio;
    sigma::Double = kCubicSplineCoefficients[dim] / h^dim;
    return CubicSplineKernel(h, dim, radius_ratio, influence_radius, sigma);
end

function kernelValue(
    r::Double, 
    kernel::CubicSplineKernel)::Double
    q::Double = r / kernel.h_;
    if q < 1.0
        return (1.0/4.0)*kernel.sigma_*(3*q*q*(q - 2) + 4);
    elseif q < 2.0
        return -1.0/4.0*kernel.sigma_*(q-2)^3;
    else
        return 0.0;
    end
end

function kernelGradient(
    r::Double, 
    kernel::CubicSplineKernel)::Double
    q::Double = r / kernel.h_;
    if q < 1.0
        return (3.0/4.0)*q*kernel.sigma_*(3*q - 4) / kernel.h_;
    elseif q < 2.0
        return -3.0/4.0*kernel.sigma_*(q-2)^2 / kernel.h_;
    else
        return 0.0;
    end
end