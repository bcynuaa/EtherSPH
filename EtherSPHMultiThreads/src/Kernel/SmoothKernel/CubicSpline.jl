#=
  @ author: bcynuaa
  @ date: 2024-01-20 23:06:33
  @ description:
 =#

struct CubicSpline{IntType <: Integer, RealType <: AbstractFloat} <: SmoothKernel
    h_::RealType
    dim_::IntType
    radius_ratio_::RealType
    influence_radius_::RealType
    sigma_::RealType
    kernel_0_::RealType
end

function kernelValue(
    r::RealType,
    kernel::CubicSpline
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
    kernel::CubicSpline
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