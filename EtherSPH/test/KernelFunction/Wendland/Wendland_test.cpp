/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-27 00:44:06
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: KernelFunction: test for Gaussian
 */

#include "KernelFunction/KernelFunction.hpp"

#include <iostream>

using namespace KernelFunction::Wendland;

int main(int argc, char const *argv[]) {
    Container::DoubleArray x_vec(2);
    const double h = 1.5e-2;
    const int dim = 2;
    const int N = 11;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Kernel Values" << std::endl;
    Container::DoubleArray kernel_value(N);
    Container::DoubleArray x_dir1 = Container::linspace<double>(0, 2*h, N);
    for (int i = 0; i < N; ++i) {
        x_vec[0] = x_dir1[i];
        kernel_value[i] = kernelValue(x_vec, h, dim);
    }
    std::cout << kernel_value << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Kernel Gradient" << std::endl;
    for (int i = 0; i < N; ++i) {
        x_vec[0] = x_dir1[i];
        std::cout << i << ": "<< kernelValueGradient(x_vec, h, dim) << std::endl;
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Kernel Second Derivatives" << std::endl;
    Container::DoubleArray kernel_laplacian(N);
    for (int i = 0; i < N; ++i) {
        x_vec[0] = x_dir1[i];
        kernel_laplacian[i] = kernelValueLaplacian(x_vec, h, dim);
    }
    std::cout << kernel_laplacian << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    return 0;
}