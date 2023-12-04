#=
  @ author: bcynuaa
  @ date: 2023-11-23 01:37:43
  @ description:
 =#

const kWendlandC4RadiusRatio::Double = 2.0;
const kWendlandC4Coefficients::AbstractVector = 
[
    5. / 8.,
    9. / 4. / pi,
    495. / 256. / pi,
];

struct WendlandC4Kernel <: AbstractKernel
    h_::Double;
    dim_::Isize;
    radius_ratio_::Double;
    influence_radius_::Double;
    sigma_::Double;
end

function WendlandC4Kernel(influence_radius::Double, dim::Isize)::WendlandC4Kernel
    radius_ratio::Double = kWendlandC4RadiusRatio;
    h::Double = influence_radius / radius_ratio;
    sigma::Double = kWendlandC4Coefficients[dim] / h^dim;
    return WendlandC4Kernel(h, dim, radius_ratio, influence_radius, sigma);
end

function kernelValue(
    r::Double, 
    kernel::WendlandC4Kernel)::Double
    q::Double = r / kernel.h_;
    if q < 2.0
        return kernel.sigma_ * (2-q)^6 * (35*q*q + 36*q + 12.) / 768.;
    else
        return 0.0;
    end
end

function kernelGradient(
    r::Double, 
    kernel::WendlandC4Kernel)::Double
    q::Double = r / kernel.h_;
    if q < 2.0
        return kernel.sigma_ * ((35. / 96.)*q*q + (7. / 48.)*q) * (q-2)^5 / kernel.h_;
    else
        return 0.0;
    end
end