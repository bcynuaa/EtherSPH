#=
  @ author: bcynuaa
  @ date: 2023-12-04 19:40:54
  @ description:
 =#

struct WendlandC2Kernel{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::WendlandC2Kernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 2.
        return kernel.sigma_ * (2. - q)^4 * (1. + 2. * q) / 16.;
    else
        return 0.;
    end
end

function kernelGradient(
    r::RealType,
    kernel::WendlandC2Kernel
)::RealType where RealType <: AbstractFloat
    q::RealType = r / kernel.h_;
    if q < 2.
        return -kernel.sigma_ / kernel.h_ * 5. / 8. * q * (2. - q)^3;
    else
        return 0.;
    end
end