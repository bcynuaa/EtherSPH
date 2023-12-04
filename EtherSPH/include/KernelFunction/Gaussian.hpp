/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 16:55:52
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: KernelFunction: head file for Gaussian
 */

#ifndef GAUSSIAN_HPP_
#define GAUSSIAN_HPP_

#include <cmath>
#include "ConstVariable/ConstVariable.hpp"
#include "Container/Container.hpp"

namespace KernelFunction {

namespace Gaussian {

const double radius_ratio = 3.0;
const Container::DoubleArray gaussian_coefficient = {1./std::sqrt(ConstVariable::Math::pi), 1./ConstVariable::Math::pi, 1./std::pow(ConstVariable::Math::pi, 3./2.)};

const double sigmaGaussian(const double h, const int dim);

const double kernelValue(const double r, const double h, const int dim);
const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff1D(const double r, const double h, const int dim);
const Container::DoubleArray kernelValueGradient(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff2D(const double r, const double h, const int dim);
const double kernelValueLaplacian(const Container::DoubleArray& x_vec, const double h, const int dim);

}; // namespace Gaussian

}; // namespace KernelFunction

#endif // GAUSSIAN_HPP_