#=
  @ author: bcynuaa
  @ date: 2023-12-04 17:28:38
  @ description:
 =#

const kCubicSplineRadiusRatio = 2.0;
const kCubicSplineCoefficients = 
[
    2.0 / 3.0,
    10. / 7. / pi,
    1. / pi,
];

struct CubicSplineKernel{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::CubicSplineKernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 1.
        return kernel.sigma_ / 4. * (3. * q^2 * (q-2.) + 4.);
    elseif q < 2.
        return kernel.sigma_ / 4. * (2. - q)^3;
    else
        return 0.;
    end
end

function kernelGradient(
    r::RealType,
    kernel::CubicSplineKernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 1.
        return kernel.sigma_ / kernel.h_ * 3. / 4. * q * (3. * q - 4.);
    elseif q < 2.
        return -kernel.sigma_ / kernel.h_ * 3. / 4. * (2. - q)^2;
    else
        return 0.;
    end
end