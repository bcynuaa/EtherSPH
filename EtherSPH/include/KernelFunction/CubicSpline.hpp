/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-26 23:47:59
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: KernelFunction: header file for CubicSpline
 */

#ifndef CUBIC_SPLINE_HPP_
#define CUBIC_SPLINE_HPP_

#include <cmath>
#include "ConstVariable/ConstVariable.hpp"
#include "Container/Container.hpp"

namespace KernelFunction {

namespace CubicSpline {

const double radius_ratio = 2.0;
const Container::DoubleArray cubic_spline_coefficient = {2./3., 10./7./ConstVariable::Math::pi, 1./ConstVariable::Math::pi};

const double sigma3(const double h, const int dim);

const double kernelValue(const double r, const double h, const int dim);
const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff1D(const double r, const double h, const int dim);
const Container::DoubleArray kernelValueGradient(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff2D(const double r, const double h, const int dim);
const double kernelValueLaplacian(const Container::DoubleArray& x_vec, const double h, const int dim);

}; // namespace CubicSpline

}; // namespace KernelFunction

#endif // CUBIC_SPLINE_HPP_