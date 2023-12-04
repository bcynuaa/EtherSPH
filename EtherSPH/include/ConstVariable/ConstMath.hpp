/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-01 21:03:48
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: ConstVariable: header file for Const Math
 */

#ifndef CONST_MATH_HPP_
#define CONST_MATH_HPP_

#include <cmath>

namespace ConstVariable {

namespace Math {

const double pi = std::acos(-1.0);

const double epsilon = 1.0e-12;

};

};

#endif // CONST_MATH_HPP_