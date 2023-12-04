/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 20:09:59
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for Wendland
 */

#ifndef WENDLAND_HPP_
#define WENDLAND_HPP_

#include <cmath>
#include "ConstVariable/ConstVariable.hpp"
#include "Container/Container.hpp"

namespace KernelFunction {

namespace Wendland {

const double radius_ratio = 2.0;
const Container::DoubleArray wendland_coefficient = {0., 7./4./ConstVariable::Math::pi, 21./16./ConstVariable::Math::pi};

const double alphaD(const double h, const int dim);

const double kernelValue(const double r, const double h, const int dim);
const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff1D(const double r, const double h, const int dim);
const Container::DoubleArray kernelValueGradient(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff2D(const double r, const double h, const int dim);
const double kernelValueLaplacian(const Container::DoubleArray& x_vec, const double h, const int dim);

};

}; // namespace KernelFunction

#endif // WENDLAND_HPP_