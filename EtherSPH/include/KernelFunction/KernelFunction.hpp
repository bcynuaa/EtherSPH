/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-23 23:47:37
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: KernelFunction: header file
 */

#ifndef KERNEL_FUNCTION_HPP_
#define KERNEL_FUNCTION_HPP_

/*
reference:
    1. https://pysph.readthedocs.io/en/latest/_modules/pysph/base/kernels.html#CubicSpline.dwdq
*/

#include "KernelFunction/CubicSpline.hpp"
#include "KernelFunction/Gaussian.hpp"
#include "KernelFunction/Wendland.hpp"
#include "KernelFunction/WendlandC4.hpp"

#endif // KERNEL_FUNCTION_HPP_