/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 16:57:53
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: source file for Gaussian
 */

#ifndef GAUSSIAN_CPP_
#define GAUSSIAN_CPP_

#include "KernelFunction/Gaussian.hpp"

namespace KernelFunction {

namespace Gaussian {

const double sigmaGaussian(const double h, const int dim) {
    return gaussian_coefficient[dim-1] / std::pow(h, dim);
};

const double kernelValue(const double r, const double h, const int dim) {
    const double q = r / h;
    if (q > 3.) {
        return 0.;
    }
    else {
        const double sigma_g = sigmaGaussian(h, dim);
        return sigma_g * std::exp(-q*q);
    }
};

const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim) {
    double r = Container::norm(x_vec);
    return kernelValue(r, h, dim);
};

const double kernelValueDiff1D(const double r, const double h, const int dim) {
    const double q = r / h;
    if (q > 3.) {
        return 0.;
    }
    else {
        const double sigma_g = sigmaGaussian(h, dim);
        return -2.0*q*sigma_g*std::exp(-q*q);
    }
};

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
    if (q > 3.) {
        return 0.;
    }
    else {
        const double sigma_g = sigmaGaussian(h, dim);
        return 2*sigma_g*(2*q*q - 1)*std::exp(-q*q);
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
        return (diff_2d/h + (dim-1)/r * diff_1d) / h;
    }
    return 0.;
};

}; // namespace Gaussian

}; // namespace KernelFunction

#endif // GAUSSIAN_CPP_