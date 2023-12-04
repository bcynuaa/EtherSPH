#=
  @ author: bcynuaa
  @ date: 2023-11-23 01:29:59
  @ description:
 =#

const kWendlandC2RadiusRatio::Double = 2.0;
const kWendlandC2Coefficients::AbstractVector = 
[
    0.,
    7. / (4. * pi),
    21. / (16. * pi),
];

struct WendlandC2Kernel <: AbstractKernel
    h_::Double;
    dim_::Isize;
    radius_ratio_::Double;
    influence_radius_::Double;
    sigma_::Double;
end

function WendlandC2Kernel(influence_radius::Double, dim::Isize)::WendlandC2Kernel
    radius_ratio::Double = kWendlandC2RadiusRatio;
    h::Double = influence_radius / radius_ratio;
    sigma::Double = kWendlandC2Coefficients[dim] / h^dim;
    return WendlandC2Kernel(h, dim, radius_ratio, influence_radius, sigma);
end

function kernelValue(
    r::Double, 
    kernel::WendlandC2Kernel)::Double
    q::Double = r / kernel.h_;
    if q < 2.0
        return kernel.sigma_ * (q-2)^4 * (2*q + 1) / 16.;
    else
        return 0.0;
    end
end

function kernelGradient(
    r::Double, 
    kernel::WendlandC2Kernel)::Double
    q::Double = r / kernel.h_;
    if q < 2.0
        return (5.0/8.0) * kernel.sigma_ * q * (q-2)^3 / kernel.h_;
    else
        return 0.0;
    end
end