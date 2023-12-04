/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-06 18:54:11
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: WendlandC4: source file
 */

#ifndef WENDLAND_C4_CPP_
#define WENDLAND_C4_CPP_

#include "KernelFunction/WendlandC4.hpp"

namespace KernelFunction {

namespace WendlandC4 {

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
        return alpha_d*pow(2 - q, 6)*(35*q*q + 36*q + 12)/768.;
    }
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
        return alpha_d*((35.0/96.0)*q*q*pow(q - 2, 5) + (7.0/48.0)*q*pow(q - 2, 5));
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
    if (q > 2.) {
        return 0.;
    }
    else {
        const double alpha_d = alphaD(h, dim);
        return alpha_d*((1.0/384.0)*alpha_d*pow(2 - q, 4)*(525*q*q + 540*q + 35*(q-2)*(q-2) + 12*(q - 2)*(35*q + 18) + 180));
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

}; // namespace WendlandC4

}; // namespace KernelFunction

#endif // WENDLAND_C4_CPP_