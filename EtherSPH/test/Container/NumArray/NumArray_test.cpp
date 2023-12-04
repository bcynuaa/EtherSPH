/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-10-20 19:58:51
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: Container: test for NumArray
 */

#include <iostream>

#include "Container/Container.hpp"

int main(int argc, char* argv[]) {
#ifdef DEBUG_MODE_
    std::cout << "DEBUG MODE" << std::endl;
#endif
    // < operator test, >, <=, >= is the same, so we only test < operator
    std::cout << "------------------------------------------------------------" << std::endl;
    Container::DoubleArray a = {2, 3, 1, -4, 7, 9, 0, 1, 2, 3};
    Container::DoubleArray b = {1, 2, 3, 4, 5, 6, 7, 8, 9, -10};
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "a < b: " << (a < b) << std::endl;
    std::cout << "a < 10.: " << (a < 10.) << std::endl;
    std::cout << "1 < a: " << (1 < a) << std::endl;
    // + and add test
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "NumArray: + and add test" << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "a + b: " << (a + b) << std::endl;
    std::cout << "a + 1: " << (a + 1) << std::endl;
    std::cout << "1 + a: " << (1 + a) << std::endl;
    a += 1.;
    std::cout << "a += 1.: " << a << std::endl;
    a += b;
    std::cout << "a += b: " << a << std::endl;
    Container::addWith(a, b);
    std::cout << "addWith(a, b): " << a << std::endl;
    Container::addWith(a, 1.);
    std::cout << "addWith(a, 1.): " << a << std::endl;
    Container::addBy(a, b);
    std::cout << "addBy(a, b): " << a << std::endl;
    Container::addBy(a, 1.);
    std::cout << "addBy(a, 1.): " << a << std::endl;

    // - and subtract test
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "NumArray: - and subtract test" << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "a - b: " << (a - b) << std::endl;
    std::cout << "a - 1: " << (a - 1) << std::endl;
    std::cout << "1 - a: " << (1 - a) << std::endl;
    a -= 1.;
    std::cout << "a -= 1.: " << a << std::endl;
    a -= b;
    std::cout << "a -= b: " << a << std::endl;
    Container::DoubleArray minus_a = -a;
    std::cout << "-a: " << minus_a << std::endl;
    Container::subtractWith(a, b);
    std::cout << "subtractWith(a, b): " << a << std::endl;
    Container::subtractWith(a, 1.);
    std::cout << "subtractWith(a, 1.): " << a << std::endl;
    Container::subtractBy(a, b);
    std::cout << "subtractBy(a, b): " << a << std::endl;
    Container::subtractBy(a, 1.);
    std::cout << "subtractBy(a, 1.): " << a << std::endl;

    // * and multiply test
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "NumArray: * and multiply test" << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "a * b: " << (a * b) << std::endl;
    std::cout << "a * 2: " << (a * 2) << std::endl;
    std::cout << "2 * a: " << (2 * a) << std::endl;
    a *= 2.;
    std::cout << "a *= 2.: " << a << std::endl;
    a *= b;
    std::cout << "a *= b: " << a << std::endl;
    Container::multiplyWith(a, b);
    std::cout << "multiplyWith(a, b): " << a << std::endl;
    Container::multiplyWith(a, 2.);
    std::cout << "multiplyWith(a, 2.): " << a << std::endl;
    Container::multiplyBy(a, b);
    std::cout << "multiplyBy(a, b): " << a << std::endl;
    Container::multiplyBy(a, 2.);
    std::cout << "multiplyBy(a, 2.): " << a << std::endl;

    // / and divide test
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "NumArray: / and divide test" << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "a / b: " << (a / b) << std::endl;
    std::cout << "a / 2: " << (a / 2) << std::endl;
    std::cout << "2 / a: " << (2 / a) << std::endl;
    a /= 2.;
    std::cout << "a /= 2.: " << a << std::endl;
    a /= b;
    std::cout << "a /= b: " << a << std::endl;
    Container::divideWith(a, b);
    std::cout << "divideWith(a, b): " << a << std::endl;
    Container::divideWith(a, 2.);
    std::cout << "divideWith(a, 2.): " << a << std::endl;
    Container::divideBy(a, b);
    std::cout << "divideBy(a, b): " << a << std::endl;
    Container::divideBy(a, 2.);
    std::cout << "divideBy(a, 2.): " << a << std::endl;

    // function library test
    std::cout << "------------------------------------------------------------" << std::endl;
    Container::DoubleArray vec(3);
    vec[0] = 3;
    vec[1] = 1;
    vec[2] = 2;
    std::cout << "vec: " << vec << std::endl;
    std::cout << "sum(vec): " << Container::sum(vec) << std::endl;
    std::cout << "prod(vec): " << Container::prod(vec) << std::endl;
    std::cout << "mean(vec): " << Container::mean(vec) << std::endl;
    std::cout << "minimum(vec): " << Container::minimum(vec) << std::endl;
    std::cout << "maximum(vec): " << Container::maximum(vec) << std::endl;
    std::cout << "linspace(0., 1., 11): " << Container::linspace(0., 1., 11) << std::endl;
    std::cout << "arange(0., 10., 1.): " << Container::arange(0., 10., 1.) << std::endl;
    std::cout << "vec . vec: " << Container::dot(vec, vec) << std::endl;
    std::cout << "norm(vec): " << Container::norm(vec) << std::endl;
    return 0;
};