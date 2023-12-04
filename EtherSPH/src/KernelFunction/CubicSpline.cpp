/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-26 23:43:40
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: KernelFunction: source file for CubicSpline
 */

#ifndef CUBIC_SPLINE_CPP_
#define CUBIC_SPLINE_CPP_

#include "KernelFunction/CubicSpline.hpp"

namespace KernelFunction {

namespace CubicSpline {

const double sigma3(const double h, const int dim) {
    return cubic_spline_coefficient[dim-1] / std::pow(h, dim);
};

const double kernelValue(const double r, const double h, const int dim) {
    const double q = r / h;
    const double sigma_3 = sigma3(h, dim);
    if (q > 2.) {
      return 0;
   }
   else if (q > 1) {
      return -1.0/4.0*sigma_3*std::pow(q - 2, 3);
   }
   else {
      return (1.0/4.0)*sigma_3*(3*q*q*(q - 2) + 4);
   }
   return 0.;
};

const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim) {
    double r = Container::norm(x_vec);
    return kernelValue(r, h, dim);
};

const double kernelValueDiff1D(const double r, const double h, const int dim) {
    const double q = r / h;
    const double sigma_3 = sigma3(h, dim);
    if (q > 2.) {
      return 0.;
   }
   else if (q > 1.) {
      return -3.0/4.0*sigma_3*std::pow(q - 2, 2);
   }
   else {
      return (3.0/4.0)*q*sigma_3*(3*q - 4);
   }
   return 0.;
}

const Container::DoubleArray kernelValueGradient(const Container::DoubleArray& x_vec, const double h, const int dim) {
    const double r = Container::norm(x_vec);
    if (r < ConstVariable::Math::epsilon) {
        return Container::zeros<double>(dim);
    }
    else {
        return (kernelValueDiff1D(r, h, dim) / r / h) * x_vec;
    }
    return Container::zeros<double>(dim);
};

const double kernelValueDiff2D(const double r, const double h, const int dim) {
    const double q = r / h;
    if (q > 2.) {
      return 0.;
   }
   else if (q > 1.) {
      const double sigma_3 = sigma3(h, dim);
      return -3.0/4.0*sigma_3*(q-2)*(q-2);
   }
   else {
      const double sigma_3 = sigma3(h, dim);
      return (3.0/2.0)*sigma_3*(3*q - 2);
   }
   return 0.;
};

const double kernelValueLaplacian(const Container::DoubleArray& x_vec, const double h, const int dim) {
    const double r = Container::norm(x_vec);
    if (r < ConstVariable::Math::epsilon) {
        return dim/h/h * kernelValueDiff2D(r, h, dim);
    }
    else {
        const double diff_1d = kernelValueDiff1D(r, h, dim);
        const double diff_2d = kernelValueDiff2D(r, h, dim);
        return (diff_1d*(dim-1)/r + diff_2d/h) / h;
    }
    return 0.;
};

}; // namespace CubicSpline

}; // namespace KernelFunction

#endif // CUBIC_SPLINE_CPP_