#=
  @ author: bcynuaa
  @ date: 2024-01-21 11:57:52
  @ description:
 =#

struct Gaussian{IntType<:Integer,RealType<:AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::Gaussian
)::RealType where {RealType<:AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 3.0
        return kernel.sigma_ * exp(-q^2)
    else
        return 0.0
    end
end

function kernelGradient(
    r::RealType,
    kernel::Gaussian
)::RealType where {RealType<:AbstractFloat}
    q::RealType = r / kernel.h_
    if q < 3.0
        return -2.0 * kernel.sigma_ / kernel.h_ * q * exp(-q^2)
    else
        return 0.0
    end
end