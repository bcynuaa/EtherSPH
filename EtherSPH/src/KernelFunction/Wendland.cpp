/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 20:09:10
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: source file for Wendland
 */

#ifndef WENDLAND_CPP_
#define WENDLAND_CPP_

#include "KernelFunction/Wendland.hpp"

namespace KernelFunction {

namespace Wendland {

const double alphaD(const double h, const int dim) {
    return wendland_coefficient[dim-1] / std::pow(h, dim);
};

const double kernelValue(const double r, const double h, const int dim) {
    const double q = r / h;
    if (q > 2.) {
        return 0.;
    }
    else {
        const double alpha_d = alphaD(h, dim);
        return alpha_d*pow(q - 2, 4)*(2*q + 1)/16.;
    }
    return 0.;
};

const double kernelValue(const Container::DoubleArray& x_vec, const double h, const int dim) {
    double r = Container::norm(x_vec);
    return kernelValue(r, h, dim);
};

const double kernelValueDiff1D(const double r, const double h, const int dim) {
    const double q = r / h;
    if (q > 2.) {
        return 0.;
    }
    else {
        const double alpha_d = alphaD(h, dim);
        return (5.0/8.0)*alpha_d*q*pow(q - 2, 3);
    }
    return 0.;
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
    if (q > 2.) {
        return 0.;
    }
    else {
        const double alpha_d = alphaD(h, dim);
        return alpha_d*(2-q)*(2-q)*(10*q - 5)/4.;
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

}; // namespace Wendland

}; // namespace KernelFunction

#endif // WENDLAND_CPP_