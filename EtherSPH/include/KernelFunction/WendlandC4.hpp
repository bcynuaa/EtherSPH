/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 18:52:43
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: WendlandC4: head file
 */

#ifndef WENDLAND_C4_HPP_
#define WENDLAND_C4_HPP_

#include <cmath>
#include "ConstVariable/ConstVariable.hpp"
#include "Container/Container.hpp"

namespace KernelFunction {

namespace WendlandC4 {

const double radius_ratio = 2.0;
const Container::DoubleArray wendland_coefficient = {5./8., 9./4./ConstVariable::Math::pi, 495./256./ConstVariable::Math::pi};

const double alphaD(const double h, const int dim);

const double kernelValue(const double r, const double h, const int dim);
const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff1D(const double r, const double h, const int dim);
const Container::DoubleArray kernelValueGradient(const Container::DoubleArray& x_vec, const double h, const int dim);

const double kernelValueDiff2D(const double r, const double h, const int dim);
const double kernelValueLaplacian(const Container::DoubleArray& x_vec, const double h, const int dim);

}; // namespace WendlandC4

}; // namespace KernelFunction

#endif // WENDLAND_C4_HPP_