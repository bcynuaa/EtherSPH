#=
  @ author: bcynuaa
  @ date: 2024-01-21 12:05:23
  @ description:
 =#

struct WendlandC4{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::WendlandC4
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 2.
        return kernel.sigma_ * (2. - q)^6 * (35. * q^2 + 36. * q + 12.) / 768.;
    else
        return 0.;
    end
end

function kernelGradient(
    r::RealType,
    kernel::WendlandC4
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 2.
        return -kernel.sigma_ / kernel.h_ * (35. / 96. * q^2 + 7. / 48. * q) * (2. - q)^5;
    else
        return 0.;
    end
end